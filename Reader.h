// Reader.h
// Loads in files
#pragma once

#include "DataStructures.h"
#include <iostream>
#include <fstream>

using namespace std;

class Reader
{
public:
    Reader();
    ~Reader();

    // Input: Name of file to be read
    // Output: If read suceeded
    bool ReadFile(string file);

    std::vector<ScalarFieldPoint> data;
    Vertex min, max;
};