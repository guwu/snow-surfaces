#include "DataStructures.h"

Vertex Vertex::operator+(const Vertex& that)
{
    Vertex result;
    result.x = this->x + that.x;
    result.y = this->y + that.y;
    result.z = this->z + that.z;
    return result;
}

Vertex Vertex::operator-(const Vertex& that)
{
    Vertex result;
    result.x = this->x - that.x;
    result.y = this->y - that.y;
    result.z = this->z - that.z;
    return result;
}

Vertex Vertex::operator-()
{
    Vertex result;
    result.x = this->x * -1;
    result.y = this->y * -1;
    result.z = this->z * -1;
    return result;
}


float Vertex::dot(Vertex& that)
{
    float d = (this->x * that.x) + (this->y * that.y) + (this->z * that.z);
    return d;
}

Vertex Vertex::cross(Vertex& that)
{
    Vertex result;
    result.x = (this->y * that.z) - (this->z * that.y);
    result.y = (this->z * that.x) - (this->x * that.z);
    result.z = (this->x * that.y) - (this->y * that.x);
    return result;
}
