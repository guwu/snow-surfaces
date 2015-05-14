// DataStructures.h
//
#ifndef DATASTRUCT
#define DATASTRUCT

#include <string>
#include <vector>
using namespace std;
class Vertex
{
public:
    Vertex() { x = 0; y = 0; z = 0; }
    Vertex operator+(const Vertex&);
    Vertex operator-(const Vertex&);
    Vertex operator-();

    float dot(Vertex&);
    Vertex cross(Vertex&);

    float x, y, z;
    float nx, ny, nz;
};

struct ScalarFieldPoint
{
    float x, y, z;
    float nx, ny, nz;
    float s;
    vector<ScalarFieldPoint*> neighbors;
};

struct Triangle
{
    Vertex *verts[3];
};

struct Edge
{
    ScalarFieldPoint *verts[2];
    bool intersection;
    Vertex *i_vert;
};

#endif