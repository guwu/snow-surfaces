// Main.cpp - Zander Clucas 2015
// 
#include <omp.h>
#include "Reader.h"
#include "DataStructures.h"
#include "MarchingCubes.h"
#include <algorithm>
#include <iostream>
#include <fstream>

// Function Prototypes
void GenScalarField(Vertex min, Vertex max);
float Kernel(float input);

// Globals
float cell_size = 0.5;
float Neighborhood = 0.5;
float part_rad = 0.4;

int numx, numy, numz;
ScalarFieldPoint ***scalar_field;
vector<Vertex> data;

// The entry point.
int main(int argc, char **argv)
{
    omp_set_num_threads(1);

    // Read File
    Reader reader;
    if (reader.ReadFile("test.txt"))
    {
        data = reader.data;
    }
    else
    {
        return 1;
    }

    // Generate Scalar Field
    GenScalarField(reader.min, reader.max);

    // Create Mesh
    MarchingCubes mc(2.f);
    vector<Triangle*> triangles = mc.triangles;
    // Output Mesh
    ofstream output;
    output.open("testmesh.out");
    if (output.is_open())
    {
        for (int i = 0; i < triangles.size(); i++)
        {
            output << "triangle { ";
            for (int j = 0; j < 3; j++)
            {
                Vertex normal, p1, p2;
                p1 = triangles[i]->verts[1] - triangles[i]->verts[0];
                p2 = triangles[i]->verts[2] - triangles[i]->verts[0];
                normal = p1.cross(p2);
                Vertex v;
                v.x = normal.x;
                v.y = normal.y;
                v.z = normal.z;

                Vertex vt;
                vt.x = triangles[i]->verts[j].x;
                vt.y = triangles[i]->verts[j].y;
                vt.z = triangles[i]->verts[j].z;
                if (vt.x < 0.0001 && vt.x > -0.0001)
                {
                    vt.x = 0.0;
                }
                if (vt.y < 0.0001 && vt.y > -0.0001)
                {
                    vt.y = 0.0;
                }
                if (vt.z < 0.0001 && vt.z > -0.0001)
                {
                    vt.z = 0.0;
                }
                output << "<" << vt.x << ", "
                    << vt.y << ", "
                    << vt.z << ">";
                if (j < 2)
                    output << ", ";
            }

            output << " }\n";
        }
        output.close();
    }
    
    return 0;
}

void GenScalarField(Vertex min, Vertex max)
{
    min.x -= 2.;
    min.y -= 2.;
    min.z -= 2.;
    max.x += 2.;
    max.y += 2.;
    max.z += 2.;

    numx = (max.x - min.x) / cell_size;
    numy = (max.y - min.y) / cell_size;
    numz = (max.z - min.z) / cell_size;

    // Allocate the scalar field
    scalar_field = new ScalarFieldPoint**[numx];
    for (int i = 0; i < numx; i++)
    {
        scalar_field[i] = new ScalarFieldPoint*[numy];
        float sx = i*cell_size;
        for (int j = 0; j < numy; j++)
        {
            float sy = j*cell_size;
            scalar_field[i][j] = new ScalarFieldPoint[numz];
            for (int k = 0; k < numz; k++)
            {
                float sz = k*cell_size;
                scalar_field[i][j][k].x = sx;
                scalar_field[i][j][k].y = sy;
                scalar_field[i][j][k].z = sz;
                scalar_field[i][j][k].s = 0;
            }
        }
    }

    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < numx; i++)
    {
        for (int j = 0; j < numy; j++)
        {
            for (int k = 0; k < numz; k++)
            {
                float s = 0;
                float rbar = 0;
                Vertex xbar;
                float ker,weight;
                float x = scalar_field[i][j][k].x;
                float y = scalar_field[i][j][k].y;
                float z = scalar_field[i][j][k].z;
                float xdist, ydist, zdist;

                for (int n = 0; n < data.size(); n++)
                {
                    Vertex node = data[n];
                    xdist = x - node.x;
                    ydist = y - node.y;
                    zdist = z - node.z;

                    if (xdist <= Neighborhood && xdist >= -Neighborhood ||
                        ydist <= Neighborhood && ydist >= -Neighborhood ||
                        zdist <= Neighborhood && zdist >= -Neighborhood)
                    {
                        ker = Kernel((1. / ((xdist * xdist) + (ydist * ydist) + (zdist * zdist))) / Neighborhood);
                        s += ker;
                    }
                    else
                    {
                        continue;
                    }
                }

                for (int n = 0; n < data.size(); n++)
                {
                    Vertex node = data[n];
                    xdist = x - node.x;
                    ydist = y - node.y;
                    zdist = z - node.z;

                    if (xdist <= Neighborhood && xdist >= -Neighborhood ||
                        ydist <= Neighborhood && ydist >= -Neighborhood ||
                        zdist <= Neighborhood && zdist >= -Neighborhood)
                    {
                        ker = Kernel((1. / ((xdist * xdist) + (ydist * ydist) + (zdist * zdist))) / Neighborhood);
                        weight = ker / s;
                        rbar += part_rad*weight;
                        xbar.x += node.x * weight;
                        xbar.y += node.y * weight;
                        xbar.z += node.z * weight;
                    }
                    else
                    {
                        continue;
                    }
                }


                xdist = x - xbar.x;
                ydist = y - xbar.y;
                zdist = z - xbar.z;

                scalar_field[i][j][k].s = (1. / ((xdist * xdist) + (ydist * ydist) + (zdist * zdist))) - rbar;
            }
        }
    }
}

float Kernel(float input)
{
    // max( 0 , (1 - s^2) ^ 3 )
    return max(0.f,pow(1-pow(input,2),3));
}

