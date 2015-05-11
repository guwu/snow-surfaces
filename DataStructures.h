// DataStructures.h
//
#ifndef DATASTRUCT
#define DATASTRUCT

#include <string>
#include <vector>

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
};

struct ScalarFieldPoint
{
    float x, y, z;
    float s;
};

struct Triangle
{
    Vertex verts[3];
};

struct Edge
{
    ScalarFieldPoint *verts[2];
    bool intersection;
    Vertex *i_vert;
};

#endif