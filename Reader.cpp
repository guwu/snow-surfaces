// Reader.cpp

#include "Reader.h"

Reader::Reader()
{

}

Reader::~Reader()
{

}

bool Reader::ReadFile(string filename)
{
    ifstream file;
    file.open(filename);
    if (file.is_open())
    {

        return true;
    }
    return false;
}