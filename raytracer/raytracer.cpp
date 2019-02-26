#include <iostream>
#include <glm/glm.hpp>
#include <SDL.h>
#include "SDLauxiliary.h"
#include "TestModelH.h"
#include <stdint.h>
#include <limits>
#include <math.h>
#include <omp.h>

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

struct Camera {
    float focalLength;
    vec4  position;
    mat4  rotation;
};

struct Light {
    vec4  position;
    vec3  color;
    vec3  indirect;
    float move;
};

#define SCREEN_WIDTH 2560
#define SCREEN_HEIGHT 1600
#define FULLSCREEN_MODE true
#define SHADOW_BIAS 0.0001f

vector<Triangle> triangles;
struct Camera camera;
struct Light  light;

bool Update();
void Draw(screen* screen);
bool TriangleIntersection(
      vec4 start,
      vec4 dir,
      Triangle triangle,
      Intersection& intersection );
bool ClosestIntersection(
      vec4 start,
      vec4 dir,
      const vector<Triangle>& triangles,
      Intersection& closestIntersection );
void RotateX( mat4& rotation, float rad );
void RotateY( mat4& rotation, float rad );
vec3 DirectLight( const Intersection& intersection );
vec3 triangleNormal(Triangle t);

int main( int argc, char* argv[] ) {
    screen *screen = InitializeSDL( SCREEN_WIDTH, SCREEN_HEIGHT, FULLSCREEN_MODE );

    // Initialise camera
    camera.focalLength = SCREEN_HEIGHT / 2;
    camera.position = vec4( 0.0, 0.0, -2.2, 1.0);
    camera.rotation = mat4(vec4(1.0, 0.0, 0.0, 1.0),
                           vec4(0.0, 1.0, 0.0, 1.0),
                           vec4(0.0, 0.0, 1.0, 1.0),
                           vec4(0.0, 0.0, 0.0, 1.0));

    light.position = vec4( 0.0f, -0.5f, -0.7f, 1.0f );
    light.color = 14.0f * vec3( 1.0, 1.0, 1.0 );
    light.indirect = 0.2f*vec3( 1, 1, 1 );
    light.move = 0.1f;

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
    #pragma omp parallel for
    for (int row = 0; row < SCREEN_HEIGHT; row++) {
        for (int col = 0; col < SCREEN_WIDTH; col++) {
            vec4 primaryRay = camera.rotation * vec4(col - SCREEN_WIDTH  / 2.0f,
                                                     row - SCREEN_HEIGHT / 2.0f,
                                                     camera.focalLength, 0.0);
            Intersection primaryIntersect;
            vec3 intersectColor = vec3(0,0,0);
            if (ClosestIntersection(camera.position, primaryRay, triangles, primaryIntersect)) {
                intersectColor = triangles[primaryIntersect.triangleIndex].color
                               * (DirectLight(primaryIntersect) + light.indirect);
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
    cout << "Render time: " << dt << " ms. " << '\t' << 1000/dt << " fps." << endl;

    SDL_Event e;
    while(SDL_PollEvent(&e)) {
        if (e.type == SDL_QUIT) {
	        return false;
	    } else if (e.type == SDL_KEYDOWN) {
	        int key_code = e.key.keysym.sym;
	        switch(key_code) {

                // Translate camera forward
                case SDLK_i:
                    camera.position += camera.rotation * vec4(0, 0, 0.1f, 0);
                    break;
                // Translate camera backwards
                case SDLK_k:
                    camera.position -= camera.rotation * vec4(0, 0, 0.1f, 0);
                    break;
                // Translate camera left
                case SDLK_j:
                    camera.position -= camera.rotation * vec4(0.1f, 0, 0, 0);
                    break;
                // Translate camera right
                case SDLK_l:
                    camera.position += camera.rotation * vec4(0.1f, 0, 0, 0);
                    break;
                // Translate camera up
                case SDLK_u:
                    camera.position -= camera.rotation * vec4(0, 0.1f, 0, 0);
                    break;
                // Translate camera down
                case SDLK_o:
                    camera.position += camera.rotation * vec4(0, 0.1f, 0, 0);
                    break;

                // Rotate camera left
	            case SDLK_LEFT:
                    RotateY(camera.rotation, -0.1f);
		            break;
                // Rotate camera right
	            case SDLK_RIGHT:
                    RotateY(camera.rotation, 0.1f);
		            break;
                // Rotate camera up
                case SDLK_UP:
                    RotateX(camera.rotation, 0.1f);
		            break;
                // Rotate camera down
	            case SDLK_DOWN:
                    RotateX(camera.rotation, -0.1f);
		            break;

                // Translate light forward
                case SDLK_w:
                    light.position += light.move * vec4(0, 0, 1.0f, 0);
                    break;
                // Translate light backwards
                case SDLK_s:
                    light.position -= light.move * vec4(0, 0, 1.0f, 0);
                    break;
                // Translate light left
                case SDLK_a:
                    light.position -= light.move * vec4(1.0f, 0, 0, 0);
                    break;
                // Translate light right
                case SDLK_d:
                    light.position += light.move * vec4(1.0f, 0, 0, 0);
                    break;
                // Translate light up
                case SDLK_q:
                    light.position -= light.move * vec4(0, 1.0f, 0, 0);
                    break;
                // Translate light down
                case SDLK_e:
                    light.position += light.move * vec4(0, 1.0f, 0, 0);
                    break;

                // Reset camera
                case SDLK_SPACE:
                    camera.position = vec4( 0.0, 0.0, -2.2, 1.0);
                    camera.focalLength = SCREEN_HEIGHT / 2;
                    camera.rotation = mat4(vec4(1.0, 0.0, 0.0, 1.0),
                                           vec4(0.0, 1.0, 0.0, 1.0),
                                           vec4(0.0, 0.0, 1.0, 1.0),
                                           vec4(0.0, 0.0, 0.0, 1.0));
                    light.position = vec4( 0.0f, -0.5f, -0.7f, 1.0f );
		            break;

                // Move camera quit
	            case SDLK_ESCAPE:
		            return false;
                    break;
            }
	    }
    }
    return true;
}

bool TriangleIntersection( vec4 start, vec4 dir, Triangle triangle, Intersection& intersection ) {

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
        intersection.distance = length(start - intersection.position);
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
    float minDistance = std::numeric_limits<float>::max();
    Intersection localIntersection;

    for( int i = 0; i < triangles.size(); i++ ) {
        if (TriangleIntersection( start, dir, triangles[i], localIntersection )) {
            found = true;
            localIntersection.triangleIndex = i;
            if (localIntersection.distance < minDistance) {
                closestIntersection = localIntersection;
                minDistance = localIntersection.distance;
            }
        }
    }
    return found;
}

void RotateX( mat4& rotation, float rad ) {
    vec4 c1(1,        0,  0,        0);
    vec4 c2(0,  cos(rad), sin(rad),      0);
    vec4 c3(0, -sin(rad), cos(rad), 0);
    vec4 c4(0,        0,  0,        1);

    mat4 R = mat4(c1, c2, c3, c4);
    rotation = R * rotation;
}

void RotateY( mat4& rotation, float rad ) {
    vec4 c1(cos(rad), 0, -sin(rad), 0);
    vec4 c2(0,        1,  0,        0);
    vec4 c3(sin(rad), 0,  cos(rad), 0);
    vec4 c4(0,        0,  0,        1);

    mat4 R = mat4(c1, c2, c3, c4);
    rotation = R * rotation;
}

vec3 DirectLight( const Intersection& intersection ) {

    vec4  normalH = triangles[intersection.triangleIndex].normal;
    vec3  normal = (vec3) normalH;
    vec3  lightVector  = light.position - intersection.position;
    float radius = length( lightVector );
    lightVector = normalize( lightVector );
    vec4 lightVectorH = vec4(lightVector.x, lightVector.y, lightVector.z, 1.0f);

    vec3 returnVector = light.color * max(dot( lightVector, normal ), 0.0f)
                      / ( (float) (4 * M_PI) * radius * radius );
    Intersection shadowIntersect;
    if (ClosestIntersection(intersection.position + normalH * SHADOW_BIAS,
                            lightVectorH, triangles, shadowIntersect)) {
        if (radius > shadowIntersect.distance) {
            returnVector = vec3(0.0f, 0.0f, 0.0f);
        }
    }
    return returnVector;
}
