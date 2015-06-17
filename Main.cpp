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
void OutputScalarField();
float* smallest(float, float*);

#ifndef CELLSIZE
#define CELLSIZE 0.4
#endif

#ifndef NEIGHBOR
#define NEIGHBOR 0.5
#endif

#ifndef RADIUS
#define RADIUS 0.1
#endif

#ifndef SURFACE
#define SURFACE .2
#endif


// Globals
float cell_size = CELLSIZE;
float Neighborhood = NEIGHBOR;
float part_rad = RADIUS;
float surface = SURFACE;
const int num_neighbors = 3;

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
    OutputScalarField();

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
    output.open("2116.pov");
    if (output.is_open())
    {
        output << "// grid cell size = " << cell_size << endl;
        output << "// Neighborhood size = " << Neighborhood << endl;
        output << "// particle radius = " << part_rad << endl;
        output << "// surface = " << surface << endl;

        output << "\n// -w1440 -h1440 +a0.1 +q11 +J0.5 +R3\n\n";
        output << "\n#version 3.7;\n\n";

        output << "global_settings{\n";
        output << "    assumed_gamma 1.0\n";
        output << "    number_of_waves 3\n";
        output << "    max_trace_level 10\n";
        output << "    //subsurface {} \n";
        output << "    //radiosity{}\n";
        output << "}\n";
        output << "\n";
        output << "#include \"colors.inc\"\n";
        output << "#include \"shapes.inc\"\n";
        output << "#include \"textures.inc\"\n";
        output << "\n";
        output << "camera{\n";
        output << "    location <58.0, 2.0, -15.0>\n";
        output << "    angle 90 //  direction z\n";
        output << "    up y\n";
        output << "    right x*image_width / image_height\n";
        output << "    look_at <0.0, 0.0, 0.0>\n";
        output << "}\n";
        output << "\n";
        output << "light_source{ <20.0, 100.0, -170.0> colour rgb <0.5, 0.5, 0.5> }\n";
        output << "light_source{ <50.0, -10.0, 10.0> colour rgb <0.1, 0.1, 0.1> }\n";
        output << "//light_source { <-35.0, 230.0, -150.0> colour White }\n";
        output << "light_source{ <90, 220, -30>\n";
        output << "    color rgb <0.5, 0.5, 0.5> * 0.1\n";
        output << "    area_light 200, 200, 10, 10\n";
        output << "    jitter\n";
        output << "}\n";
        output << "\n";
        output << "#declare Snow =\n";
        output << "texture{\n";
        output << "    pigment{ color rgb<.38, .75, 1> }\n";
        output << "    finish{\n";
        output << "        ambient 0.1\n";
        output << "        diffuse 0.7\n";
        output << "        brilliance 0.8\n";
        output << "        specular 0.4\n";
        output << "        roughness 0.04\n";
        output << "    }\n";
        output << "    normal{\n";
        output << "        granite 0.2\n";
        output << "        scale 0.3\n";
        output << "    }\n";
        output << "    finish{\n";
        output << "        brilliance 0.75\n";
        output << "        phong 0.2\n";
        output << "        phong_size 5\n";
        output << "        subsurface{ translucency <0.3, 0.41, 0.48> }\n";
        output << "\n";
        output << "        //emission .2\n";
        output << "        //use with  radiosity instead\n";
        output << "    }\n";
        output << "}\n";
        output << "\n";
        output << "#declare Asphalt = texture{\n";
        output << "pigment{ color rgb<1, 1, 1>*0.2 }\n";
        output << "normal{ bumps 0.5 scale 0.005 }\n";
        output << "finish{ diffuse 0.9 phong 0.1 }\n";
        output << "}\n\n";

        output << "fog{\n";
        output << "    distance 800 color rgb <0.39, 0.51, 0.61>\n";
        output << "}\n";
        output << "\n";
        output << "background{ color rgb <0.39, 0.51, 0.61> }\n";
        output << "\n"; 
        output << "/* Ground plane */\n";
        output << "plane{\n";
        output << "    y, -60\n";
        output << "    texture{\n";
        output << "    Asphalt\n";
        output << "}\n";
        output << "}\n\n";

        output << "mesh {\n";
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

        output << "\ntexture {\n\tSnow\n}interior { ior 1.31 }\n}";

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


void OutputScalarField()
{
    ofstream output;
    std::string sfile = "scalarfield." + std::to_string(cell_size) + "." + std::to_string(Neighborhood) + "." + std::to_string(num_neighbors) + "." + std::to_string(part_rad) + ".field";
    output.open(sfile, fstream::out);
    if (output.is_open())
    {
        output << "// grid cell size = " << cell_size << endl;
        output << "// Neighborhood size = " << Neighborhood << endl;
        output << "// num neighbors = " << num_neighbors << endl;
        output << "// particle radius = " << part_rad << endl;
        
        output << "// " << numx << " " << numy << " " << numz << " " << endl << endl;
        for (int i = 0; i < numx; i++)
        {
            for (int j = 0; j < numy; j++)
            {
                for (int k = 0; k < numz; k++)
                {
                    output << i << " " << j << " " << k << " ";
                    output << " " << scalar_field[i][j][k].x;
                    output << " " << scalar_field[i][j][k].y;
                    output << " " << scalar_field[i][j][k].z;
                    output << " " << scalar_field[i][j][k].s << endl;
                }
            }
        }
    }
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
    
    // Create Neighborhoods
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
                    #pragma omp critical
                    scalar_field[i][j][k].neighbors.push_back(&data[n]);
                }
            }
        }
        
    }

    // Generate Scalar Field
    ScalarFieldPoint *neighbors_array[num_neighbors];
    float neighbors_dist_array[num_neighbors];
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
                
                for (int m = 0; m < num_neighbors; m++)
                {
                    neighbors_dist_array[m] = 9999999.f;
                }
                
                int neighbors = 0;

                // minimize the neighborhood
                for (int n = 0; n < scalar_field[i][j][k].neighbors.size(); n++)
                {
                    ScalarFieldPoint *node = scalar_field[i][j][k].neighbors[n];
                    xdist = x - node->x;
                    ydist = y - node->y;
                    zdist = z - node->z;
                    dist = (sqrt((xdist * xdist) + (ydist * ydist) + (zdist * zdist)));
                    
                    if (dist <= Neighborhood)
                    {
                        // Cheap implementation, just take the first n neighbors.
                        neighbors_array[neighbors] = node;
                        neighbors_dist_array[neighbors] = dist;
                        neighbors++;
                        if (neighbors == num_neighbors)
                        {
                            break;
                        }
                    }
                    else
                    {
                        continue;
                    }
                }

                // Determine the kernel sum denominator
                for (int n = 0; n < neighbors; n++)
                {
                    ScalarFieldPoint *node = neighbors_array[n];
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
                 
                // Determine the kernel, xbar and rbar
                for (int n = 0; n < neighbors; n++)
                {
                    ScalarFieldPoint *node = neighbors_array[n];
                    xdist = x - node->x;
                    ydist = y - node->y;
                    zdist = z - node->z;
                    dist = (sqrt((xdist * xdist) + (ydist * ydist) + (zdist * zdist)));

                    if (dist <= Neighborhood)
                    {
                        ker = Kernel(dist / Neighborhood);
                        weight = (ker / s)*(rand()/RAND_MAX)*100.f; // random bumpyness
                        rbar += part_rad*weight*(rand() / RAND_MAX)*5.f;
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

                scalar_field[i][j][k].s = (sqrt((xdist * xdist) + (ydist * ydist) + (zdist * zdist))) - rbar;
            }
        }
    }


    // Generate Normals
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
    
    return max((double)0., (double)pow(1 - pow(input, 2), 3));
    //if(input>=0)
    //    return max((double)0., (double)(.03*sin(40 * input)*input + 2 - 2.022*pow(input, 2)) / 2.);
	//else
    //    return max((double)0., (double)(cos(20 * input)*.1 + 1.9 - 2.153*pow(input, 2) + 0.159*fabs((double)input)) / 2.);
}

