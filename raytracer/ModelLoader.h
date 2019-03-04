#ifndef MODEL_LOADER
#define MODEL_LOADER

#define LINES 1024

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <glm/glm.hpp>
#include <SDL.h>
#include <limits.h>
#include "TestModelH.h"

using namespace std;
using glm::vec3;
using glm::mat3;
using glm::vec4;
using glm::mat4;

bool LoadModel( std::vector<Triangle>& triangles, const char* filename ) {
    cout << "al doilea segfault" << endl;
    char* folder = (char*) malloc(8 * sizeof(char));
    char  path[256];
    folder  = "Models/";
    path[0] = '\0';
    strcat(path, folder);
    strcat(path, filename);
    ifstream objfile;
    objfile.open(path);
    if( objfile.fail() ) {
        cerr << "Error: " << strerror(errno) << endl;
        return false;
    }

    string line;
    string header, s[3];
    int n = 0, vInd = 0, vtInd = 0, vnInd = 0;
    vec3 vertices[LINES];
    vec3 textures[LINES];
    vec3  normals[LINES];
    while( !objfile.eof() ) {
        getline (objfile, line);
        istringstream input( line );
        input >> header >> s[0] >> s[1] >> s[2];

        if        (header == "v") {
            input >> header >> s[0] >> s[1] >> s[2];
            vertices[vInd].x = strtof((s[0]).c_str(), 0);
            vertices[vInd].y = strtof((s[1]).c_str(), 0);
            vertices[vInd].z = strtof((s[2]).c_str(), 0);
            cout << "v " << vertices[vInd].x << "\t" << vertices[vInd].y << "\t" << vertices[vInd].z << endl;
            vInd++;
        } else if (header == "vt") {
            // to implement
        } else if (header == "vn") {
            input >> header >> s[0] >> s[1] >> s[2];
            normals[vnInd].x = strtof((s[0]).c_str(), 0);
            normals[vnInd].y = strtof((s[1]).c_str(), 0);
            normals[vnInd].z = strtof((s[2]).c_str(), 0);
            cout << "vn " << normals[vnInd].x << "\t" << normals[vnInd].y << "\t" << normals[vnInd].z << endl;
            vnInd++;
        } else if (header == "f") {
            vec3 color = vec3(2.0f, 2.0f, 2.0f);
            vec4 v[3], normals[3];
            int vertexInd[3], textureInd[3], normalInd[3];
            input >> header >> s[0] >> s[1] >> s[2];

            for ( int pos = 0; pos < 3; pos++) {
                /* v1/vt1/vn1  v2/vt2/vn2  v3/vt3/vn3
                ** 0 -> vertices "v"
                ** 1 -> textures "vt"
                ** 2 -> normals  "vn"
                */
                vector<string> fields;
                string::const_iterator end     = s[pos].end();
                string::const_iterator current = s[pos].begin();
                string::const_iterator next = find( current, end, '/' );
                while ( next != end ) {
                    fields.push_back( string( current, next ) );
                    current = next + 1;
                    next = find( current, end, '/' );
                }
                fields.push_back( string( current, next ) );

                vertexInd[pos]  = (fields[0] == "")   ? -1 : stoi((fields[0]).c_str());
                textureInd[pos] = (fields[1] == "" ||
                                   fields.size() < 2) ? -1 : stoi((fields[1]).c_str());
                normalInd[pos]  = (fields[2] == ""||
                                   fields.size() < 3) ? -1 : stoi((fields[2]).c_str());

                vec4 vertex = vec4( vertices[vertexInd[pos] - 1].x,
                                    vertices[vertexInd[pos] - 1].y,
                                    vertices[vertexInd[pos] - 1].z,
                                    1.0f);
                v[pos] = vertex;
            }

            Triangle triangle = Triangle( v[0], v[1], v[2], color );
            triangles.push_back(triangle);
        } else {
            n++;
        }
    }
    objfile.close();
    return true;
}
#endif
