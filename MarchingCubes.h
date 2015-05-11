// MarchingCubes.h
//
#pragma once

#include "DataStructures.h"
#include <omp.h>

using namespace std;

extern float cell_size;
extern float Neighborhood;
extern float part_rad;
extern int numx, numy, numz;
extern ScalarFieldPoint ***scalar_field;

class MarchingCubes
{
public:
    MarchingCubes(float st);
    vector<Triangle*> triangles;

private:
    float Surface;
    void MarchingCube(int x, int y, int z);

};

