#include <string>
#include <iostream>
#include <fstream>
#include <vector>

#include <maya/MPoint.h>
#include <maya/MVector.h>

#include "WireSculptPlugin.h"

using namespace std;

// checks to make sure file is obj 
bool WireSculptPlugin::GetFileExtension(const string& filePath) {
    size_t dotPosition = filePath.find_last_of(".");
    return (filePath.substr(dotPosition + 1) == "obj");
}

// split string helper 
vector<string> WireSculptPlugin::SplitString(const string& input, char delimiter) {
    vector<string> tokens;
    size_t start = 0;
    size_t end = input.find(delimiter);

    while (end != string::npos) {
        tokens.push_back(input.substr(start, end - start));
        start = end + 1;
        end = input.find(delimiter, start);
    }

    tokens.push_back(input.substr(start));

    return tokens;
}

// processes the obj file 
// returns false: issue with file (contents or type) 
// returns true: success 
bool WireSculptPlugin::ProcessFile(std::string filePath) {
    // incorrect obj type 
    if (!GetFileExtension(filePath)) {
        return false; 
    }

    // failure to open 
    ofstream fin;
    fin.open(filePath);
    if (!fin.is_open()) {
        return false;
    }

    string line;

    // temp storage 
    vector<MVector> pos; 
    vector<MVector> vertexNormals;

    while (getline(cin, line)) {
        if (line == "-1") {
            break;
        }
        // adding vert positions to temp vector
        else if (line.substr(0) == "v") {
            vector<string> fields = SplitString(line, ' ');
            pos.push_back(MVector(stof(fields[1]), stof(fields[2]), stof(fields[3])));
        }

        // calculating normals by going through faces 
        else if (line.substr(0) == "f") {

            vector<string> fields = SplitString(line, ' ');
            // temp storage to calc normals w/ vert pos
            vector<float> faceVerts; 

            // going through verticies that make up a face 
            for (int i = 1; i < fields.size(); i++) {
                vector<string> triplets = SplitString(line, '/');

                faceVerts.push_back(stof(triplets[0]));
            }

            // draw edges with verts to get cross prod
            for (int i = 0; i < faceVerts.size(); i++) {
                
                MVector current = pos[i];
                MVector prev = pos[(i - 1)%pos.size()];
                MVector next = pos[(i + 1)%pos.size()];

                MVector crossProd = (prev - current) ^ (next-current);
                crossProd.normalize();

                // for each calc can make Vertex obj and add to list
                verticies.push_back(Vertex(MPoint(pos[0]), crossProd));
            }
        }
    }

    fin.close();
    return true;
}

std::vector<Vertex>* WireSculptPlugin::GetVerticies() {
    return &(this->verticies);
}