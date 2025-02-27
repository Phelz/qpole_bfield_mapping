using namespace std;

#include <iostream>
#include <omp.h>
#include <chrono>


#include "constants.h"
#include "vector_utils.h"  // Include the header where vector operations are defined FIRST
#include "biot_savart_routines.h"



void biot_savart_chukman(const PRECISION_TYPE* s0,
    const PRECISION_TYPE* s1,
    const PRECISION_TYPE* r,
    PRECISION_TYPE* B) {

    PRECISION_TYPE i[3], k0[3], k1[3], A[3];
    PRECISION_TYPE norm_i, norm_k0, norm_k1, norm_A_sqr, a0;

    // i = s1 - s0
    subtract(s1, s0, i);
    norm_i = norm(i);
    if (norm_i < TOO_SMALL) return;
    scale(i, 1.0 / norm_i, i); // Normalize i

    subtract(s0, r, k0); // k0 = s0 - r
    subtract(s1, r, k1); // k1 = s1 - r

    norm_k0 = norm(k0);
    norm_k1 = norm(k1);
    if (norm_k0 < TOO_SMALL || norm_k1 < TOO_SMALL) return;

    // A = cross(k0, i)
    cross_product(k0, i, A);
    norm_A_sqr = dot_product(A, A);
    if (norm_A_sqr < TOO_SMALL) return;

    scale(A, 1.0 / norm_A_sqr, A); 
    scale(k0, 1.0 / norm_k0, k0); 
    scale(k1, 1.0 / norm_k1, k1); 


    // a0 = dot(k1 - k0, i)
    PRECISION_TYPE k1_minus_k0[3];
    subtract(k1, k0, k1_minus_k0);
    a0 = dot_product(k1_minus_k0, i);

    // Accumulate the result: B = B + A * a0
    PRECISION_TYPE scaled_A[3];
    scale(A, a0, scaled_A); // Scale A by a0
    add(B, scaled_A);       // Add scaled A to B (accumulation)

}


void calc_bfield_parallel_chukman(const uint32_t num_seg,
    const PRECISION_TYPE* seg_start,
    const PRECISION_TYPE* seg_end,
    const uint32_t num_grid,
    const PRECISION_TYPE* grid_xyz,
    PRECISION_TYPE* bfield_xyz,
    bool print_progress) {

    // Reset magnetic field grid
    fill(bfield_xyz, bfield_xyz + 3 * num_grid, 0.0);

    auto start_time = chrono::high_resolution_clock::now();

    for (int seg_indx = 0; seg_indx < num_seg; seg_indx++) {
        const PRECISION_TYPE* start_point = &seg_start[3 * seg_indx];
        const PRECISION_TYPE* end_point = &seg_end[3 * seg_indx];


        try
        {
            #pragma omp parallel for
            for (int grid_indx = 0; grid_indx < num_grid; grid_indx++) {
                biot_savart_chukman(start_point, end_point, &grid_xyz[3 * grid_indx], &bfield_xyz[3 * grid_indx]);
            }
        }
        catch (std::exception& ex) {
            std::cout << "Simulation failed." << std::endl;
            std::cout << ex.what() << std::endl;
            return;

        }

        if (print_progress) {if (seg_indx % (num_seg / 10) == 0 || seg_indx == num_seg - 1) { // Update for every 10% of progress or at the last segment
                int percentage = static_cast<int>((static_cast<double>(seg_indx) / num_seg) * 100);
                cout << "Progress: " << percentage << "% (" << seg_indx << " / " << num_seg << ")" << endl;
            }
        }
    }
    
    auto end_time = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed = end_time - start_time;
    cout << "Bfield Calculation finished (Chukman). Time elapsed: " << elapsed.count() << " seconds" << endl;
}

