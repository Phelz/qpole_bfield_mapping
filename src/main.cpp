using namespace std;

// #include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <omp.h>
#include <sstream>
#include <map>
// #include <windows.h> // ! No longer using windows

#include "constants.h"
//#include "vector_utils.h"
#include "ParticleTrace.h"
#include "biot_savart_routines.h"


void printMemoryUsage() {
  system("wmic OS get FreePhysicalMemory,TotalVisibleMemorySize");
}

void generate_uniform_grid(vector<PRECISION_TYPE> &grid_xyz,
                           PRECISION_TYPE x_min, PRECISION_TYPE x_max,
                           PRECISION_TYPE y_min, PRECISION_TYPE y_max,
                           PRECISION_TYPE z_min, PRECISION_TYPE z_max,
                           uint32_t num_points) {

  grid_xyz.clear();
  grid_xyz.reserve(3 * num_points * num_points * num_points);

  PRECISION_TYPE x_step = (x_max - x_min) / (num_points - 1);
  PRECISION_TYPE y_step = (y_max - y_min) / (num_points - 1);
  PRECISION_TYPE z_step = (z_max - z_min) / (num_points - 1);

  for (uint32_t i = 0; i < num_points; i++) {
    PRECISION_TYPE x = x_min + i * x_step;
    for (uint32_t j = 0; j < num_points; j++) {
      PRECISION_TYPE y = y_min + j * y_step;
      for (uint32_t k = 0; k < num_points; k++) {
        PRECISION_TYPE z = z_min + k * z_step;
        grid_xyz.push_back(x);
        grid_xyz.push_back(y);
        grid_xyz.push_back(z);
      }
    }
  }
}


void save_magnetic_field(const string &filename, const PRECISION_TYPE *grid_xyz,
                         const PRECISION_TYPE *bfield_xyz,
                         uint32_t num_points) {

  // Open the file
  ofstream outfile(filename, ios::out | ios::trunc);
  if (!outfile.is_open()) {
    cerr << "Error: Could not open file " << filename << " for writing!"
         << endl;
    return;
  }

  // Set precision
  outfile << fixed << setprecision(OUTPUT_PRECISION);

  // Write the header
  outfile << "x,y,z,Bx,By,Bz\n";

  // Write data
  for (uint32_t i = 0; i < num_points; i++) {
    outfile << grid_xyz[3 * i] << "," << grid_xyz[3 * i + 1] << ","
            << grid_xyz[3 * i + 2] << "," << bfield_xyz[3 * i] << ","
            << bfield_xyz[3 * i + 1] << "," << bfield_xyz[3 * i + 2] << "\n";
  }

  outfile.close();
  cout << "Magnetic field data saved to " << filename << endl;
}

void generate_straight_wire(uint32_t num_segments, PRECISION_TYPE zmin,
                            PRECISION_TYPE dz,
                            vector<PRECISION_TYPE> &segment_start,
                            vector<PRECISION_TYPE> &segment_end) {

  segment_start.clear();
  segment_end.clear();
  segment_start.reserve(3 * num_segments);
  segment_end.reserve(3 * num_segments);

  for (uint32_t seg = 0; seg < num_segments; seg++) {
    segment_start.push_back(0.0);
    segment_start.push_back(0.0);
    segment_start.push_back(zmin + seg * dz);

    segment_end.push_back(0.0);
    segment_end.push_back(0.0);
    segment_end.push_back(zmin + (seg + 1) * dz);
  }
}

void generate_circular_loop(vector<PRECISION_TYPE> &segment_start,
                            vector<PRECISION_TYPE> &segment_end,
                            PRECISION_TYPE radius, uint32_t num_segments) {

  segment_start.clear();
  segment_end.clear();
  segment_start.reserve(3 * num_segments);
  segment_end.reserve(3 * num_segments);

  PRECISION_TYPE dtheta = 2.0 * PI / num_segments;

  for (uint32_t i = 0; i < num_segments; i++) {
    PRECISION_TYPE theta1 = i * dtheta;
    PRECISION_TYPE theta2 = (i + 1) * dtheta;

    segment_start.push_back(radius * cos(theta1));
    segment_start.push_back(radius * sin(theta1));
    segment_start.push_back(0.0);

    segment_end.push_back(radius * cos(theta2));
    segment_end.push_back(radius * sin(theta2));
    segment_end.push_back(0.0);
  }
}

void test_circular_loop(const uint32_t num_segments, PRECISION_TYPE loop_radius,
                        const uint32_t grid_size) {
  // Create Loop
  vector<PRECISION_TYPE> segment_start(num_segments * 3);
  vector<PRECISION_TYPE> segment_end(num_segments * 3);
  generate_circular_loop(segment_start, segment_end, loop_radius, num_segments);

  // Grid
  vector<PRECISION_TYPE> grid(grid_size * grid_size * grid_size * 3);
  generate_uniform_grid(grid, -loop_radius * 5, loop_radius * 5,
                        -loop_radius * 5, loop_radius * 5, -loop_radius * 5,
                        loop_radius * 5, grid_size);

  // Initialize magnetic field vectors
  vector<PRECISION_TYPE> bfield_chukman(grid.size(), 0.0);

  // Pass raw pointers using .data()
  calc_bfield_parallel_chukman(num_segments, segment_start.data(),
                               segment_end.data(), grid.size() / 3, grid.data(),
                               bfield_chukman.data(), true);

  // Save magnetic field data
  save_magnetic_field("../data/B_field_loop_chukman.txt", grid.data(),
                      bfield_chukman.data(), grid.size() / 3);
}

void test_wire_z(const uint32_t num_segments, PRECISION_TYPE zmin,
                 PRECISION_TYPE zmax, const uint32_t grid_size) {
  // Create Wire
  PRECISION_TYPE dz = (zmax - zmin) / num_segments;
  vector<PRECISION_TYPE> segment_start(num_segments * 3);
  vector<PRECISION_TYPE> segment_end(num_segments * 3);
  generate_straight_wire(num_segments, zmin, dz, segment_start, segment_end);

  // Grid
  vector<PRECISION_TYPE> grid(grid_size * grid_size * grid_size * 3);
  generate_uniform_grid(grid, ZMIN, ZMAX, ZMIN, ZMAX, ZMIN, ZMAX, grid_size);

  // Initialize magnetic field vectors
  vector<PRECISION_TYPE> bfield_chukman(grid.size(), 0.0);

  // Pass raw pointers using .data()
  calc_bfield_parallel_chukman(num_segments, segment_start.data(),
                               segment_end.data(), grid.size() / 3, grid.data(),
                               bfield_chukman.data(), true);

  // Save magnetic field data

  save_magnetic_field("../data/B_field_wire_chukman.txt", grid.data(),
                      bfield_chukman.data(), grid.size() / 3);
}

void test_particle_trace(const uint32_t grid_size) {


  // Create grid
  vector<PRECISION_TYPE> grid;
  generate_uniform_grid(grid, XMIN - X_RANGE / 2, XMAX + X_RANGE / 2,
                        	  YMIN - Y_RANGE / 2, YMAX + Y_RANGE / 2,
                        	  ZMIN - Z_RANGE / 2, ZMAX + Z_RANGE / 2, grid_size);

  string directory = TRACES_DIR;
  string filepath = directory + "/particle_1.csv";
  ParticleTrace particle_trace(filepath);

  // Retrieve segments
  const auto &segment_start = particle_trace.getSegmentStart();
  const auto &segment_end = particle_trace.getSegmentEnd();
  uint32_t num_segments = particle_trace.getNumSegments();

  // Initialize magnetic field vectors
  vector<PRECISION_TYPE> bfield_chukman(grid.size(), 0.0);

  // Pass raw pointers using .data()
  calc_bfield_parallel_chukman(num_segments, segment_start, segment_end,
                               grid.size() / 3, grid.data(),
                               bfield_chukman.data(), true);
  
  // scale the magnetic field with MU_0/4PI
  for (size_t i = 0; i < bfield_chukman.size(); i++) {
    bfield_chukman[i] *= MU0 / (4.0 * PI);
  }

  // Save magnetic field data
  save_magnetic_field("../data/B_field_wire_chukman.txt", grid.data(),
                      bfield_chukman.data(), grid.size() / 3);
}

void test_all_particle_traces(const uint32_t grid_size) {

  // Create grid
  vector<PRECISION_TYPE> grid;
  generate_uniform_grid(grid, XMIN - X_RANGE / 2, XMAX + X_RANGE / 2,
                            YMIN - Y_RANGE / 2, YMAX + Y_RANGE / 2,
                            ZMIN - Z_RANGE / 2, ZMAX + Z_RANGE / 2, grid_size);

  string directory = TRACES_DIR;

  // Initialize magnetic field vectors to accumulate results
  vector<PRECISION_TYPE> bfield_chukman(grid.size(), 0.0);

  // Read currents from the file
  map<int, PRECISION_TYPE> particle_currents;
  ifstream currents_file("../data/particle_currents.csv");
  if (currents_file.is_open()) {
    string line;
    while (getline(currents_file, line)) {
      istringstream iss(line);
      int particle_number;
      PRECISION_TYPE current;
      char comma;
      if (iss >> particle_number >> comma >> current) {
        particle_currents[particle_number] = current;
      }
    }
    currents_file.close();
  } else {
    cerr << "Error: Could not open particle currents file!" << endl;
    return;
  }

  // Iterate over all particle trace files (particle_0.csv to particle_599.csv)
  for (int i = 0; i < 600; i++) {

    cout << "Processing particle trace #: " << i << "..." << endl;
    string filepath = directory + "/particle_" + to_string(i) + ".csv";
    // First, see / the path exists
    if (!ifstream(filepath).good()) {
      cout << "Skipping missing particle trace #: " << i << endl;
      continue; // Skip this particle
    }

    // Get the current for this particle
    PRECISION_TYPE current = particle_currents[i];
    // Lets do some print statements for testing
    cout << "Current for particle " << i << " is " << current << endl;


    ParticleTrace particle_trace(filepath);

    // Retrieve segments
    const auto &segment_start = particle_trace.getSegmentStart();
    const auto &segment_end = particle_trace.getSegmentEnd();
    uint32_t num_segments = particle_trace.getNumSegments();

    // * No need anymore
    // // Check if the particle is frozen
    // if (particle_trace.checkFrozen(y_midpoint, z_midpoint, x_range)) {
    // 	cout << "Skipping frozen particle trace #: " << i << endl;
    // 	continue;  // Skip this particle
    // }

    // Pass raw pointers using .data()
    vector<PRECISION_TYPE> bfield_chukman_temp(grid.size(), 0.0);

    calc_bfield_parallel_chukman(num_segments, segment_start, segment_end,
                                 grid.size() / 3, grid.data(),
                                 bfield_chukman_temp.data(), false);


    for (size_t i = 0; i < bfield_chukman_temp.size(); i++) {
      bfield_chukman_temp[i] *= (MU0*current) / (4.0 * PI);
    }

    // Accumulate B-field results for this particle trace
    for (size_t grid_indx = 0; grid_indx < bfield_chukman.size(); grid_indx++) {
      bfield_chukman[grid_indx] += bfield_chukman_temp[grid_indx];
    }

    if (i % 30 == 0) { // 600 traces, 10% is every 60th trace
      cout << "Progress: " << (i / 6) << "% done." << endl;
      printMemoryUsage();
    }
  }

  // Save the final accumulated magnetic field data
  save_magnetic_field("B_field_traces_chukman_total.txt", grid.data(),
                      bfield_chukman.data(), grid.size() / 3);

  cout << "Finished processing all particle traces. Total B-field accumulated "
          "and saved."
       << endl;
}

int main() {
  cout << fixed << setprecision(OUTPUT_PRECISION);
  omp_set_dynamic(0);  // Disable dynamic adjustment of threads
  omp_set_num_threads(NUM_THREADS);


  cout << "Grid bounds: " << XMIN << " " << XMAX << " " << YMIN << " "
       << YMAX << " " << ZMIN << " " << ZMAX << endl;


  // SetPriorityClass(GetCurrentProcess(), HIGH_PRIORITY_CLASS);

  // test_wire_z(100, -10.0, 10.0, 100);
  // test_circular_loop(100, 2.0, 100);
  // test_particle_trace(30);
  test_all_particle_traces(30);

  cout << "This is working" << endl;
}
