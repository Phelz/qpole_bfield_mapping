
// extern "C" {
//     void calc_bfield_parallel_se(const uint32_t num_seg,
//         const PRECISION_TYPE* seg_start,
//         const PRECISION_TYPE* seg_end,
//         const uint32_t num_grid,
//         const PRECISION_TYPE* grid_xyz,
//         PRECISION_TYPE* bfield_xyz,
//         bool print_progress);
// }


// B-field Parallel for SE
void calc_bfield_parallel_se(
    const uint32_t num_seg,          // number of wire segments
    const PRECISION_TYPE *seg_start, // start of segments
    const PRECISION_TYPE *seg_end,   // end of segments
    const uint32_t num_grid,         // number of grid points
    const PRECISION_TYPE *grid_xyz,  // the grid
    PRECISION_TYPE *bfield_xyz, bool print_progress);

    
void biot_savart_se(const PRECISION_TYPE *start, const PRECISION_TYPE *end,
    const PRECISION_TYPE *grid, PRECISION_TYPE *B);