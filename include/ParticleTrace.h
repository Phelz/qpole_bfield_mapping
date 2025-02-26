#pragma once

using namespace std;

#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <algorithm>

#include "constants.h"

class ParticleTrace {

    private:
    PRECISION_TYPE* segment_start;  // Flat array for segment start points
    PRECISION_TYPE* segment_end;    // Flat array for segment end points
    PRECISION_TYPE current;
    size_t num_segments;

public:
    ParticleTrace(const string& filename);
    ~ParticleTrace();

    void readCSV(const string& filename);
    // bool checkFrozen(PRECISION_TYPE y_midpoint, PRECISION_TYPE z_midpoint, PRECISION_TYPE x_range);

    const PRECISION_TYPE* getSegmentStart() { return segment_start; };
    const PRECISION_TYPE* getSegmentEnd() { return segment_end; };
    const size_t getNumSegments() const { return num_segments; };
};

ParticleTrace::ParticleTrace(const string& filename)
    : segment_start(nullptr), segment_end(nullptr), num_segments(0) {
    readCSV(filename);
}

ParticleTrace::~ParticleTrace() {
    delete[] segment_start;
    delete[] segment_end;
}



void ParticleTrace::readCSV(const string& filename) {
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Error: Could not open file " << filename << endl;
        return;
    }

    string line;
    vector<PRECISION_TYPE> positions;

    // Read the first line (header)
    getline(file, line);

    // Read the data
    while (getline(file, line)) {
        stringstream ss(line);
        string value;
        int col_indx = 0;

        // Read x, y, z
        while (getline(ss, value, ',')) {
            if (col_indx == 3) break;
            positions.push_back(stold(value));
            col_indx++;
        }
    }
    file.close();

    // Calculate number of segments
    num_segments = positions.size() / 3 - 1;

    // Allocate memory for the segment start and end arrays
    segment_start = new PRECISION_TYPE[num_segments * 3];
    segment_end   = new PRECISION_TYPE[num_segments * 3];

    for (size_t i = 0; i < num_segments; i++) {
        segment_start[i * 3 + 0] = positions[i * 3 + 0];
        segment_start[i * 3 + 1] = positions[i * 3 + 1];
        segment_start[i * 3 + 2] = positions[i * 3 + 2];
        segment_end[i * 3 + 0]   = positions[(i + 1) * 3 + 0];
        segment_end[i * 3 + 1]   = positions[(i + 1) * 3 + 1];
        segment_end[i * 3 + 2]   = positions[(i + 1) * 3 + 2];
    }
}


