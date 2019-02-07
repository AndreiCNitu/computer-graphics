#include <iostream>
#include <glm/glm.hpp>
#include <SDL.h>
#include "SDLauxiliary.h"
#include "TestModelH.h"
#include <stdint.h>
#include <limits>
#include <math.h>

using namespace std;
using glm::vec3;
using glm::mat3;
using glm::vec4;
using glm::mat4;

struct Intersection {
    vec4 position;
    float distance;
    int triangleIndex;
};

#define SCREEN_WIDTH 1280
#define SCREEN_HEIGHT 800
#define FULLSCREEN_MODE true

SDL_Event event;
vector<Triangle> triangles;
float focalLength = SCREEN_HEIGHT / 2;
vec4 cameraPos( 0.0, 0.0, -3.0, 1.0);

bool Update();
void Draw(screen* screen);
bool TriangleIntersection( vec4 start, vec4 dir, int i, Triangle triangle, Intersection& intersection );
bool ClosestIntersection(
      vec4 start,
      vec4 dir,
      const vector<Triangle>& triangles,
      Intersection& closestIntersection );

int main( int argc, char* argv[] ) {

    screen *screen = InitializeSDL( SCREEN_WIDTH, SCREEN_HEIGHT, FULLSCREEN_MODE );

    LoadTestModel( triangles );

    while ( Update() ) {
        Draw(screen);
        SDL_Renderframe(screen);
    }

    SDL_SaveImage( screen, "screenshot.bmp" );

    KillSDL(screen);
    return 0;
}

void Draw(screen* screen) {
    /* Clear buffer */
    memset(screen->buffer, 0, screen->height * screen->width * sizeof(uint32_t));

    // Loop throught all pixels
    for (int row = 0; row < SCREEN_HEIGHT; row++) {
        for (int col = 0; col < SCREEN_WIDTH; col++) {
            vec4 ray(col - SCREEN_WIDTH / 2.0f, row - SCREEN_HEIGHT / 2.0f, focalLength, 0.0);
            Intersection intersect;
            vec3 intersectColor = vec3(0, 0, 0);
            if (ClosestIntersection(cameraPos, ray, triangles, intersect)) {
                intersectColor = triangles[intersect.triangleIndex].color;
            }
            PutPixelSDL(screen, col, row, intersectColor);
        }
    }
}

bool Update() {
    static int t = SDL_GetTicks();
    /* Compute frame time */
    int t2 = SDL_GetTicks();
    float dt = float(t2-t);
    t = t2;

    SDL_Event e;
    while(SDL_PollEvent(&e)) {
        if (e.type == SDL_QUIT) {
	        return false;
	    } else if (e.type == SDL_KEYDOWN) {
	        int key_code = e.key.keysym.sym;
	        switch(key_code) {
	            case SDLK_UP:
		            /* Move camera forward */
		            break;
	            case SDLK_DOWN:
		            /* Move camera backwards */
		            break;
	            case SDLK_LEFT:
		            /* Move camera left */
		            break;
	            case SDLK_RIGHT:
		            /* Move camera right */
		            break;
	            case SDLK_ESCAPE:
		            /* Move camera quit */
		            return false;
                    break;
            }
	    }
    }
    return true;
}

bool TriangleIntersection( vec4 start, vec4 dir, int i, Triangle triangle, Intersection& intersection ) {

    vec4 v0 = triangle.v0;
    vec4 v1 = triangle.v1;
    vec4 v2 = triangle.v2;

    // v0 + u e1 + v e2 = s + t d
    vec3 e1 = vec3(v1.x - v0.x, v1.y - v0.y, v1.z - v0.z); // e1 = v1 - v0
    vec3 e2 = vec3(v2.x - v0.x, v2.y - v0.y, v2.z - v0.z); // e2 = v2 - v0
    vec3 b  = vec3(start.x - v0.x, start.y - v0.y, start.z - v0.z); // s - v0

    // v0 + ue1 + ve2 = s + td
    mat3 A( -dir, e1, e2 );
    vec3 x  = glm::inverse( A ) * b;
    float t = x.x;
    float u = x.y;
    float v = x.z;

    intersection.distance = std::numeric_limits<float>::max();

    if (t >= 0 && u >= 0 && v >= 0 && (u + v) <= 1 ) {
        intersection.position = start + t * dir;
        float dx2 = (start.x - intersection.position.x) * (start.x - intersection.position.x);
        float dy2 = (start.y - intersection.position.y) * (start.y - intersection.position.y);
        float dz2 = (start.z - intersection.position.y) * (start.z - intersection.position.z);
        intersection.distance = sqrt(dx2 + dy2 + dz2);
        intersection.triangleIndex = i;
        return true;
    } else {
        return false;
    }
}

bool ClosestIntersection(
    vec4 start,
    vec4 dir,
    const vector<Triangle>& triangles,
    Intersection& closestIntersection ) {

    bool found = false;
    float max = std::numeric_limits<float>::max();
    Intersection localIntersection;

    for( int i = 0; i < triangles.size(); i++ ) {
        if (TriangleIntersection( start, dir, i, triangles[i], localIntersection )) {
            found = true;
            if (localIntersection.distance < max) {
                closestIntersection = localIntersection;
                max = localIntersection.distance;
            }
        }
    }
    return found;
}
