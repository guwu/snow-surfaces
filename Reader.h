// Reader.h
// Loads in files

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

    std::vector<Vertex> data;
};