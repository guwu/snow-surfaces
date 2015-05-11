// Reader.cpp

#include "Reader.h"

Reader::Reader()
{
    min.x = 9999;
    min.y = 9999;
    min.z = 9999;

    max.x = -9999;
    max.y = -9999;
    max.z = -9999;
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
        string line;
        string povout;
        string s1, s2;
        double time = 0;
        float x, y, z;
        double sxx, sxy, syy, syz, szz, sxz, mass;
        int matnum;
        

        // Header (hard coded bs)
        getline(file, line); // Source
        file >> s1 >> time >> s2; // Time
        getline(file, line); // (return char)
        getline(file, line); // Export
        getline(file, line); // Data
        getline(file, line); // Format
        getline(file, line); // End
        if (!line.compare("EndHeader") == 0)
        {
            printf("Output file lacking header, unable to determine time.\n");
            file.close();
            return false;
        }

        // Data
        while (file >> x >> y >> z
            >> sxx >> syy >> szz
            >> sxy >> syz >> sxz >> matnum >> mass)
        {
            Vertex v;
            v.x = x;
            v.y = y;
            v.z = z;
            data.push_back(v);
            if (x < min.x)
                min.x = x;
            if (x > max.x)
                max.x = x;
            if (y < min.y)
                min.y = y;
            if (y > max.y)
                max.y = y;
            if (z < min.z)
                min.z = z;
            if (z > max.z)
                max.z = z;
        }
        file.close();
        return true;
    }
    return false;
}