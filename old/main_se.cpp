
  // Read grid bounds
  PRECISION_TYPE x_min, x_max, y_min, y_max, z_min, z_max, x_range, x_midpoint,
      y_range, y_midpoint, z_range, z_midpoint;

  if (!read_grid_bounds("../data/constants_test21.csv", x_min, x_max, x_range,
                        x_midpoint, y_min, y_max, y_range, y_midpoint, z_min,
                        z_max, z_range, z_midpoint)) {
    throw runtime_error("Failed to read grid bounds from constants_test21.csv");
  }

// Function to read min/max values from CSV file
bool read_grid_bounds(const string &filename, PRECISION_TYPE &x_min,
    PRECISION_TYPE &x_max, PRECISION_TYPE &x_range,
    PRECISION_TYPE &x_midpoint, PRECISION_TYPE &y_min,
    PRECISION_TYPE &y_max, PRECISION_TYPE &y_range,
    PRECISION_TYPE &y_midpoint, PRECISION_TYPE &z_min,
    PRECISION_TYPE &z_max, PRECISION_TYPE &z_range,
    PRECISION_TYPE &z_midpoint) {

ifstream file(filename);
if (!file.is_open()) {
cerr << "Error: Could not open file " << filename << endl;

// Check if the file exists or is accessible
if (errno == ENOENT) {
cerr << "Error: The file does not exist." << endl;
} else if (errno == EACCES) {
cerr << "Error: Permission denied to open the file." << endl;
} else {
cerr << "Error: Unknown error occurred. Error code: " << errno << ""
<< endl;
}

return false;
}

string line;
getline(file, line); // Read the first line (header)
getline(file, line); // Read the data line

stringstream ss(line);
char comma; // For handling commas in CSV

// Read values in order
ss >> x_max >> comma >> y_max >> comma >> z_max >> comma >> x_min >> comma >>
y_min >> comma >> z_min >> comma >> x_range >> comma >> y_range >>
comma >> z_range >> comma >> x_midpoint >> comma >> y_midpoint >> comma >>
z_midpoint;

file.close();
return true;
}

// vector<PRECISION_TYPE> bfield_se(grid.size(), 0.0);

// calc_bfield_parallel_se(num_segments, segment_start.data(), segment_end.data(),
// grid.size() / 3, grid.data(), bfield_se.data(), true);


// save_magnetic_field("../data/B_field_loop_se.txt", grid.data(), bfield_se.data(), grid.size() / 3);

// vector<PRECISION_TYPE> bfield_se(grid.size(), 0.0);
 
// calc_bfield_parallel_se(num_segments, segment_start.data(), segment_end.data(),
// grid.size() / 3, grid.data(), bfield_se.data(), true);

// save_magnetic_field("../data/B_field_wire_se.txt", grid.data(), bfield_se.data(), grid.size() / 3);


// vector<PRECISION_TYPE> bfield_se(grid.size(), 0.0);
// calc_bfield_parallel_se(num_segments, segment_start, segment_end, grid.size() / 3, grid.data(), bfield_se.data(), true);
// save_magnetic_field("../data/B_field_wire_se.txt", grid.data(), bfield_se.data(), grid.size() / 3);

// vector<PRECISION_TYPE> bfield_se(grid.size(), 0.0);
// vector<PRECISION_TYPE> bfield_se_temp(grid.size(), 0.0);
// calc_bfield_parallel_se(num_segments, segment_start, segment_end, grid.size() / 3, grid.data(), bfield_se_temp.data(), false);
// bfield_se[grid_indx] += bfield_se_temp[grid_indx];
// save_magnetic_field("B_field_traces_se_total.txt", grid.data(), bfield_se.data(), grid.size() / 3);