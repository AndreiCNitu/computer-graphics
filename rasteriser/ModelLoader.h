#ifndef MODEL_LOADER
#define MODEL_LOADER

#define LINES 128 * 1024

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <glm/glm.hpp>
#include <SDL.h>
#include <limits.h>
#include "TestModelH.h"

using namespace std;
using glm::vec2;
using glm::vec3;
using glm::vec4;
using glm::mat4;

bool LoadModel( vector<Triangle>& triangles, const char* filename );
void RotateTriangles( vector<Triangle>& triangles );
void ScaleTriangles(  vector<Triangle>& triangles );
bool LoadCornellBox(  vector<Triangle>& triangles );

bool LoadModel( vector<Triangle>& triangles, const char* filename ) {
    char* folder = (char*) malloc(8 * sizeof(char));
    char  path[256];
    folder  = "../Models/";
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
    vec2 textures[LINES];
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
            vInd++;
        } else if (header == "vt") {
            input >> header >> s[0] >> s[1] >> s[2];
            textures[vtInd].x = strtof((s[0]).c_str(), 0);
            textures[vtInd].y = strtof((s[1]).c_str(), 0);
        } else if (header == "vn") {
            input >> header >> s[0] >> s[1] >> s[2];
            normals[vnInd].x = strtof((s[0]).c_str(), 0);
            normals[vnInd].y = strtof((s[1]).c_str(), 0);
            normals[vnInd].z = strtof((s[2]).c_str(), 0);
            vnInd++;
        } else if (header == "f") {
            // TODO: set color
            vec3 color( 0.20f, 0.20f, 0.90f );
            vec4 v[3], uv[3], normals[3];
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

                vec2 uvCoord = vec2( textures[textureInd[pos] - 1].x,
                                     textures[textureInd[pos] - 1].y);
                uv[pos] = vertex;
            }

            Triangle triangle = Triangle( v[2], v[1], v[0], color );
            // Triangle triangle = Triangle( v[2], v[1], v[0], uv[0], uv[1], uv[2] );
            triangles.push_back(triangle);
        } else {
            n++;
        }
    }
    objfile.close();

    RotateTriangles( triangles );
    ScaleTriangles(  triangles );
    return true;
}

void RotateTriangles( vector<Triangle>& triangles) {
    // Rotate around X-axis
    // (for some reason X is flipped in the obj file)
    for( size_t i = 0; i < triangles.size(); ++i ) {
        mat4 Rx = mat4(vec4(1,  0,  0, 0),
                       vec4(0, -1,  0, 0),
                       vec4(0,  0, -1, 0),
                       vec4(0,  0,  0, 1));
        triangles[i].v0 = Rx * triangles[i].v0;
    	triangles[i].v1 = Rx * triangles[i].v1;
    	triangles[i].v2 = Rx * triangles[i].v2;
    }
}

void ScaleTriangles( vector<Triangle>& triangles ) {
    float MAX  = 0.0f; // max coords
    float yMin;        // min y coord
    yMin = numeric_limits<float>::max();

    for( size_t i = 0; i < triangles.size(); ++i ) {
        if        (triangles[i].v0.x > MAX) {
            MAX =  triangles[i].v0.x;
        } else if (triangles[i].v1.x > MAX) {
            MAX =  triangles[i].v1.x;
        } else if (triangles[i].v2.x > MAX) {
            MAX =  triangles[i].v2.x;
        } else if (triangles[i].v0.y > MAX) {
            MAX =  triangles[i].v0.y;
        } else if (triangles[i].v1.y > MAX) {
            MAX =  triangles[i].v1.y;
        } else if (triangles[i].v2.y > MAX) {
            MAX =  triangles[i].v2.y;
        } else if (triangles[i].v0.z > MAX) {
            MAX =  triangles[i].v0.z;
        } else if (triangles[i].v1.z > MAX) {
            MAX =  triangles[i].v1.z;
        } else if (triangles[i].v2.z > MAX) {
            MAX =  triangles[i].v2.z;
        }

        if        (triangles[i].v0.y < yMin) {
            yMin =  triangles[i].v0.y;
        } else if (triangles[i].v1.y < yMin) {
            yMin =  triangles[i].v1.y;
        } else if (triangles[i].v2.y < yMin) {
            yMin =  triangles[i].v2.y;
        }
    }
    yMin *= 1 / (2.0f * MAX);

    for( size_t i = 0; i < triangles.size(); ++i ) {
        // Scale object to fit inside box
    	triangles[i].v0 *= 1 / (2.0f * MAX);
    	triangles[i].v1 *= 1 / (2.0f * MAX);
    	triangles[i].v2 *= 1 / (2.0f * MAX);

        // Move object towards camera
        triangles[i].v0.z -= 0.6f;
    	triangles[i].v1.z -= 0.6f;
    	triangles[i].v2.z -= 0.6f;

        triangles[i].v0.y += (0.5f - yMin);
    	triangles[i].v1.y += (0.5f - yMin);
    	triangles[i].v2.y += (0.5f - yMin);

    	triangles[i].v0.w = 1.0f;
    	triangles[i].v1.w = 1.0f;
    	triangles[i].v2.w = 1.0f;

    	triangles[i].ComputeNormal();
	}
}

bool LoadCornellBox( vector<Triangle>& triangles ) {

    // Defines colors:
    vec3 red(    0.75f, 0.15f, 0.15f );
    vec3 yellow( 0.75f, 0.75f, 0.15f );
    vec3 green(  0.15f, 0.75f, 0.15f );
    vec3 cyan(   0.15f, 0.75f, 0.75f );
    vec3 blue(   0.15f, 0.15f, 0.75f );
    vec3 purple( 0.75f, 0.15f, 0.75f );
    vec3 white(  0.75f, 0.75f, 0.75f );

    triangles.clear();
    triangles.reserve( 5*2*3 );

    // ---------------------------------------------------------------------------
    // Room

    float L = 555;			// Length of Cornell Box side.

    vec4 A(L,0,0,1);
    vec4 B(0,0,0,1);
    vec4 C(L,0,L,1);
    vec4 D(0,0,L,1);

    vec4 E(L,L,0,1);
    vec4 F(0,L,0,1);
    vec4 G(L,L,L,1);
    vec4 H(0,L,L,1);

    // Floor:
    triangles.push_back( Triangle( C, B, A, green ) );
    triangles.push_back( Triangle( C, D, B, green ) );

    // Left wall
    triangles.push_back( Triangle( A, E, C, purple ) );
    triangles.push_back( Triangle( C, E, G, purple ) );

    // Right wall
    triangles.push_back( Triangle( F, B, D, yellow ) );
    triangles.push_back( Triangle( H, F, D, yellow ) );

    // Ceiling
    triangles.push_back( Triangle( E, F, G, cyan ) );
    triangles.push_back( Triangle( F, H, G, cyan ) );

    // Back wall
    triangles.push_back( Triangle( G, D, C, white ) );
    triangles.push_back( Triangle( G, H, D, white ) );

    // ----------------------------------------------
    // Scale to the volume [-1,1]^3

	for( size_t i=0; i<triangles.size(); ++i ) {
    	triangles[i].v0 *= 2/L;
    	triangles[i].v1 *= 2/L;
    	triangles[i].v2 *= 2/L;

    	triangles[i].v0 -= vec4(1,1,1,1);
    	triangles[i].v1 -= vec4(1,1,1,1);
    	triangles[i].v2 -= vec4(1,1,1,1);

    	triangles[i].v0.w = 1.0;
    	triangles[i].v1.w = 1.0;
    	triangles[i].v2.w = 1.0;

    	triangles[i].ComputeNormal();
	}

    return true;
}

#endif
