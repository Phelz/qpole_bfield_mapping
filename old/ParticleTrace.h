
// * No need to checkFrozen anymore
// bool ParticleTrace::checkFrozen(PRECISION_TYPE y_midpoint, PRECISION_TYPE z_midpoint, PRECISION_TYPE x_range) {
//     if (num_segments == 0) return true; // No segments, consider frozen

//     // Given starting box vertices
//     PRECISION_TYPE box_vertices[4][3] = {
//         {0.1034488454, 0.0757665179, 0.03999450579},
//         {0.1034488454, 0.0745003179, 0.03999450579},
//         {0.1034488454, 0.0757665179, 0.1299945058},
//         {0.1034488454, 0.0745003179, 0.1299945058}
//     };

//     // Compute reflected end box
//     PRECISION_TYPE end_box[4][3];
//     for (int i = 0; i < 4; i++) {
//         end_box[i][0] = box_vertices[i][0];  // x remains the same
//         end_box[i][1] = -(box_vertices[i][1] - y_midpoint) + y_midpoint;  // Reflect y
// 		end_box[i][2] = -(box_vertices[i][2] - z_midpoint) + z_midpoint;  // Reflect z
//     }

//     // Compute bounding box for end region
//     PRECISION_TYPE x_tol = x_range / 5;
//     PRECISION_TYPE x_min_end = end_box[0][0] - x_tol;
//     PRECISION_TYPE x_max_end = end_box[0][0] + x_tol;

//     PRECISION_TYPE y_min_end = min(min(end_box[0][1], end_box[1][1]), min(end_box[2][1], end_box[3][1]));
//     PRECISION_TYPE y_max_end = max(max(end_box[0][1], end_box[1][1]), max(end_box[2][1], end_box[3][1]));

//     PRECISION_TYPE z_min_end = min(min(end_box[0][2], end_box[1][2]), min(end_box[2][2], end_box[3][2]));
//     PRECISION_TYPE z_max_end = max(max(end_box[0][2], end_box[1][2]), max(end_box[2][2], end_box[3][2]));


//     // Get final segment position
//     PRECISION_TYPE x_final = segment_end[(num_segments - 1) * 3 + 0];
//     PRECISION_TYPE y_final = segment_end[(num_segments - 1) * 3 + 1];
//     PRECISION_TYPE z_final = segment_end[(num_segments - 1) * 3 + 2];

//     // Check if final position is within the bounds
//     return !(x_min_end <= x_final && x_final <= x_max_end &&
//         y_min_end <= y_final && y_final <= y_max_end &&
//         z_min_end <= z_final && z_final <= z_max_end);
// }
