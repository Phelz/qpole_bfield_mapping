#pragma once

#include "constants.h"

void biot_savart_chukman(const PRECISION_TYPE *s0, const PRECISION_TYPE *s1,
                         const PRECISION_TYPE *r, PRECISION_TYPE *B);


// B-field Parallel for Chunkman
void calc_bfield_parallel_chukman(
    const uint32_t num_seg,          // number of wire segments
    const PRECISION_TYPE *seg_start, // start of segments
    const PRECISION_TYPE *seg_end,   // end of segments
    const uint32_t num_grid,         // number of grid points
    const PRECISION_TYPE *grid_xyz,  // the grid
    PRECISION_TYPE *bfield_xyz,
    bool print_progress); // Output B-field grid
