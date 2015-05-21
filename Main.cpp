// Main.cpp - Zander Clucas 2015
// 
#include <omp.h>
#include "Reader.h"
#include "DataStructures.h"
#include "MarchingCubes.h"
#include <algorithm>
#include <iostream>
#include <fstream>
#include <math.h>

// Function Prototypes
void GenScalarField(Vertex min, Vertex max);
float Kernel(float input);
void OutputPovRay();
void OutputRenderman();

// Globals
float cell_size = 0.4;
float Neighborhood = 1.25f;
float part_rad = 0.1;
float surface = 1.25;

int numx, numy, numz;
ScalarFieldPoint ***scalar_field;
vector<ScalarFieldPoint> data;
vector<Triangle*> triangles;

// The entry point.
int main(int argc, char **argv)
{
    omp_set_num_threads(1);

    // Read File
    Reader reader;
    if (reader.ReadFile("2116"))
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
    MarchingCubes mc; 
    triangles = mc.March(surface);

    // Output Mesh
    OutputPovRay();
    
    return 0;
}


void OutputPovRay()
{
    ofstream output;
    output.open("testmesh.pov");
    if (output.is_open())
    {
        for (int i = 0; i < triangles.size(); i++)
        {
            output << "smooth_triangle { ";
            for (int j = 0; j < 3; j++)
            {
                Vertex vt;
                vt.x = triangles[i]->verts[j]->x;
                vt.y = triangles[i]->verts[j]->y;
                vt.z = triangles[i]->verts[j]->z;
                vt.nx = triangles[i]->verts[j]->nx;
                vt.ny = triangles[i]->verts[j]->ny;
                vt.nz = triangles[i]->verts[j]->nz;
/*
                p1.x = triangles[i]->verts[1]->x - triangles[i]->verts[0]->x;
                p1.y = triangles[i]->verts[1]->y - triangles[i]->verts[0]->y;
                p1.z = triangles[i]->verts[1]->z - triangles[i]->verts[0]->z;
                p2.x = triangles[i]->verts[2]->x - triangles[i]->verts[0]->x;
                p2.y = triangles[i]->verts[2]->y - triangles[i]->verts[0]->y;
                p2.z = triangles[i]->verts[2]->z - triangles[i]->verts[0]->z;
                normal = p1.cross(p2);

                normal.x = triangles.*/

                
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
                output << ",<" << vt.nx << ", "
                    << vt.ny << ", "
                    << vt.nz << ">";
                if (j < 2)
                    output << ", ";
            }

            output << " }\n";
        }
        output.close();
    }
}

void OutputRenderman()
{
    //fstream blobby;
    //std::string sfile = "2116.rib";
    //blobby.open(sfile, fstream::out);
    //if (blobby.is_open())
    //{
    //    blobby << "GeneralPolygon " << (int)points.size() << " [\n";

    //    for (int i = 0; i < (int)points.size(); i++)
    //    {
    //        blobby << "\t1001 " << i * 16 << "\t\t\t\t\t# " << i << endl;
    //        //blobby << "\t1004 " << "0 0 0 0 0\t\t\t\t\t# " << i << endl;
    //    }
    //    blobby << "\t0 " << (int)points.size();
    //    for (int i = 0; i < (int)points.size(); i++)
    //    {
    //        blobby << " " << i;
    //    }
    //    blobby << "\t\t# " << (int)points.size() + 1 << endl;
    //    blobby << "] [\n";
    //    for (int i = 0; i < (int)points.size(); i++)
    //    {
    //        blobby << "\t1 0 0 0  0 1 0 0  0 0 1 0  ";
    //        blobby << points[i].x << " " << points[i].y << " " << points[i].z;
    //        blobby << " 1\t\t\t# " << i << endl;
    //    }
    //    blobby << "] [ \"${RMANTREE}/etc/impl_cube.dll\" ]\n";
    //    blobby.close();
    //}

}

void GenScalarField(Vertex min, Vertex max)
{
    min.x -= 5.;
    min.y -= 5.;
    min.z -= 5.;
    max.x += 5.;
    max.y += 5.;
    max.z += 5.;

    numx = (max.x - min.x) / cell_size;
    numy = (max.y - min.y) / cell_size;
    numz = (max.z - min.z) / cell_size;


    double time0 = omp_get_wtime();
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
                scalar_field[i][j][k].x = min.x+sx;
                scalar_field[i][j][k].y = min.y+sy;
                scalar_field[i][j][k].z = min.z+sz;
                scalar_field[i][j][k].s = 0;
            }
        }
    }
    
    #pragma omp parallel for schedule(dynamic)
    for (int n = 0; n < data.size(); n++)
    {
        ScalarFieldPoint node = data[n];
        float x = (node.x - min.x) / cell_size;
        float y = (node.y - min.y) / cell_size;
        float z = (node.z - min.z) / cell_size;

        float Neigborsteps = (Neighborhood*2.) / cell_size; // times 1.5 for a buffer

        for (int i = (int)floor(x - Neigborsteps); i < (int)ceil(x + Neigborsteps); i++)
        {
            for (int j = (int)floor(y - Neigborsteps); j < (int)ceil(y + Neigborsteps); j++)
            {
                for (int k = (int)floor(z - Neigborsteps); k < (int)ceil(z + Neigborsteps); k++)
                {
                    scalar_field[i][j][k].neighbors.push_back(&data[n]);
                }
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
                ScalarFieldPoint xbar;
                xbar.x = 0;
                xbar.y = 0;
                xbar.z = 0;
                float ker,weight;
                float x = scalar_field[i][j][k].x;
                float y = scalar_field[i][j][k].y;
                float z = scalar_field[i][j][k].z;
                float xdist, ydist, zdist, dist;

                for (int n = 0; n < scalar_field[i][j][k].neighbors.size(); n++)
                {
                    ScalarFieldPoint *node = scalar_field[i][j][k].neighbors[n];
                    xdist = x - node->x;
                    ydist = y - node->y;
                    zdist = z - node->z;
                    dist = (sqrt((xdist * xdist) + (ydist * ydist) + (zdist * zdist)));


                    if (dist <= Neighborhood)
                    {
                        ker = Kernel( dist / Neighborhood);
                        s += ker;
                    }
                    else
                    {
                        continue;
                    }
                }
                if (s == 0)
                {
                    scalar_field[i][j][k].s = 0;
                    continue;
                }

                for (int n = 0; n < scalar_field[i][j][k].neighbors.size(); n++)
                {
                    ScalarFieldPoint *node = scalar_field[i][j][k].neighbors[n];
                    xdist = x - node->x;
                    ydist = y - node->y;
                    zdist = z - node->z;
                    dist = (sqrt((xdist * xdist) + (ydist * ydist) + (zdist * zdist)));

                    if (dist <= Neighborhood)
                    {
                        ker = Kernel(dist / Neighborhood);
                        weight = ker / s;
                        rbar += dist;// part_rad*weight;
                        xbar.x += node->x * weight;
                        xbar.y += node->y * weight;
                        xbar.z += node->z * weight;
                    }
                    else
                    {
                        continue;
                    }
                }

                xdist = x - xbar.x;
                ydist = y - xbar.y;
                zdist = z - xbar.z;

                scalar_field[i][j][k].s = (sqrt((xdist * xdist) + (ydist * ydist) + (zdist * zdist))) - part_rad;
            }
        }
    }


    // Normals
#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < numx; i++)
    {
        for (int j = 0; j < numy; j++)
        {
            for (int k = 0; k < numz; k++)
            {
                //nx
                if (i == 0)
                {
                    scalar_field[i][j][k].nx = -(scalar_field[i][j][k].s - scalar_field[i + 1][j][k].s) / (2 * cell_size);
                }
                else if (i == numx-1)
                {
                    scalar_field[i][j][k].nx = -(scalar_field[i - 1][j][k].s - scalar_field[i][j][k].s) / (2 * cell_size);
                }
                else
                {
                    scalar_field[i][j][k].nx = -(scalar_field[i - 1][j][k].s - scalar_field[i + 1][j][k].s) / (2 * cell_size);
                }

                //ny
                if (j == 0)
                {
                    scalar_field[i][j][k].ny = -(scalar_field[i][j][k].s - scalar_field[i][j + 1][k].s) / (2 * cell_size);
                }
                else if (j == numy - 1)
                {
                    scalar_field[i][j][k].ny = -(scalar_field[i][j - 1][k].s - scalar_field[i][j][k].s) / (2 * cell_size);
                }
                else
                {
                    scalar_field[i][j][k].ny = -(scalar_field[i][j - 1][k].s - scalar_field[i][j + 1][k].s) / (2 * cell_size);
                }

                //nz
                if (k == 0)
                {
                    scalar_field[i][j][k].nz = -(scalar_field[i][j][k].s - scalar_field[i][j][k + 1].s) / (2 * cell_size);
                }
                else if (k == numz - 1)
                {
                    scalar_field[i][j][k].nz = -(scalar_field[i][j][k - 1].s - scalar_field[i][j][k].s) / (2 * cell_size);
                }
                else
                {
                    scalar_field[i][j][k].nz = -(scalar_field[i][j][k - 1].s - scalar_field[i][j][k + 1].s) / (2 * cell_size);
                }
            }
        }
    }

    double time1 = omp_get_wtime();
    cout << (time1 - time0) << endl;
}

float Kernel(float input)
{
	// max( 0 , (1 - s^2) ^ 3 )
    
	return max(0.,pow(1-pow(input,2),3));
//if(input>=0)
	//   return max(0.,(.03*sin(40*input)*input+2-2.022*pow(input,2))/2);
	//else
	//   return max(0.,(cos(20*input)*.1+1.9-2.153*pow(input,2)+0.159*abs(input))/2);
}

