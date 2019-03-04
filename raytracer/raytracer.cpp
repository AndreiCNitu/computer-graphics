#include <iostream>
#include <glm/glm.hpp>
#include <SDL.h>
#include "SDLauxiliary.h"
#include "TestModelH.h"
#include "ModelLoader.h"
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

struct Camera {
    float focalLength;
    vec4  position;
    mat4  rotation;
    float moveSpeed;
};

struct Light {
    vec4  position;
    vec3  color;
    vec3  indirect;
    float moveSpeed;
};

#define SCREEN_WIDTH  2560
#define SCREEN_HEIGHT 1600
#define FULLSCREEN_MODE true
#define SHADOW_BIAS 0.0001f
/*  6   2   7
**    \ | /
** 3 -- 1 -- 4
**    / | \
**  8   5   9
** Set to 1/5/9 */
#define AA_SAMPLES 1


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
void InitialiseParams(); // Initialise camera and light

int main( int argc, char* argv[] ) {

    if (argc > 3) {
        cout << "Too many arguments." << endl;
        cout << "Usage: ./raytracer --cornell-box or ./raytracer --load <model.obj>" << endl;
        exit(1);
    } else if (argc == 2 && ( strcmp("--cornell-box", argv[1]) == 0 ) ||
               argc == 1) {
        LoadTestModel( triangles );
        cout << "Cornell box loaded successfully." << endl;
    } else if (argc == 3 && strcmp("--load",        argv[1]) == 0 ) {
        if (LoadModel( triangles, argv[2] )) {
            cout << argv[2] << " model loaded successfully." << endl;
        } else {
            cout << "Unable to load " << argv[2] << endl;
            exit(1);
        }
    } else {
        cout << "Unknown command." << endl;
        cout << "Usage: ./raytracer --cornell-box or ./raytracer --load <model.obj>" << endl;
        exit(1);
    }

    screen *screen = InitializeSDL( SCREEN_WIDTH, SCREEN_HEIGHT, FULLSCREEN_MODE );
    InitialiseParams();

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
    float col_offset[9] = {0.0f,  0.0f, -0.5f, 0.5f, 0.0f, -0.5f,  0.5f, -0.5f, 0.5f};
    float row_offset[9] = {0.0f, -0.5f,  0.0f, 0.0f, 0.5f, -0.5f, -0.5f,  0.5f, 0.5f};

    // Loop throught all pixels
    #pragma omp parallel for
    for (int row = 0; row < SCREEN_HEIGHT; row++) {
        for (int col = 0; col < SCREEN_WIDTH; col++) {
            vec3 pixelColor = vec3(0,0,0);
            for(int sample = 0; sample < AA_SAMPLES; sample++) {
                vec4 primaryRay = camera.rotation * vec4(col - SCREEN_WIDTH  / 2.0f + col_offset[sample],
                                                         row - SCREEN_HEIGHT / 2.0f + row_offset[sample],
                                                         camera.focalLength, 0.0);
                Intersection primaryIntersect;
                if (ClosestIntersection(camera.position, primaryRay, triangles, primaryIntersect)) {
                    pixelColor += triangles[primaryIntersect.triangleIndex].color
                                * (DirectLight(primaryIntersect) + light.indirect);
                }
            }
            PutPixelSDL(screen, col, row, pixelColor / (float) AA_SAMPLES);
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
                    camera.position += camera.moveSpeed * camera.rotation * vec4(0, 0, 0.1f, 0);
                    break;
                // Translate camera backwards
                case SDLK_k:
                    camera.position -= camera.moveSpeed * camera.rotation * vec4(0, 0, 0.1f, 0);
                    break;
                // Translate camera left
                case SDLK_j:
                    camera.position -= camera.moveSpeed * camera.rotation * vec4(0.1f, 0, 0, 0);
                    break;
                // Translate camera right
                case SDLK_l:
                    camera.position += camera.moveSpeed * camera.rotation * vec4(0.1f, 0, 0, 0);
                    break;
                // Translate camera up
                case SDLK_u:
                    camera.position -= camera.moveSpeed * camera.rotation * vec4(0, 0.1f, 0, 0);
                    break;
                // Translate camera down
                case SDLK_o:
                    camera.position += camera.moveSpeed * camera.rotation * vec4(0, 0.1f, 0, 0);
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
                    light.position += light.moveSpeed * vec4(0, 0, 1.0f, 0);
                    break;
                // Translate light backwards
                case SDLK_s:
                    light.position -= light.moveSpeed * vec4(0, 0, 1.0f, 0);
                    break;
                // Translate light left
                case SDLK_a:
                    light.position -= light.moveSpeed * vec4(1.0f, 0, 0, 0);
                    break;
                // Translate light right
                case SDLK_d:
                    light.position += light.moveSpeed * vec4(1.0f, 0, 0, 0);
                    break;
                // Translate light up
                case SDLK_q:
                    light.position -= light.moveSpeed * vec4(0, 1.0f, 0, 0);
                    break;
                // Translate light down
                case SDLK_e:
                    light.position += light.moveSpeed * vec4(0, 1.0f, 0, 0);
                    break;

                // Reset camera and light
                case SDLK_SPACE:
                    InitialiseParams();
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

void InitialiseParams() {
    camera.focalLength = SCREEN_HEIGHT;
    camera.position = vec4( 0.0, 0.0, -3.001f, 1.0);
    camera.rotation = mat4(vec4(1.0, 0.0, 0.0, 1.0),
                           vec4(0.0, 1.0, 0.0, 1.0),
                           vec4(0.0, 0.0, 1.0, 1.0),
                           vec4(0.0, 0.0, 0.0, 1.0));
    camera.moveSpeed = 1.0f;

    light.position = vec4( 0.0f, -0.5f, -0.7f, 1.0f );
    light.color = 14.0f * vec3( 1.0, 1.0, 1.0 );
    light.indirect = 0.2f*vec3( 1, 1, 1 );
    light.moveSpeed = 0.1f;
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
    const float detA  = glm::determinant( A );

    float t = glm::determinant( mat3(   b, e1, e2 ) ) / detA;
    if (t < 0) {
        return false;
    }
    float u = glm::determinant( mat3( -dir, b, e2 ) ) / detA;
    if (u < 0) {
        return false;
    }
    float v = glm::determinant( mat3( -dir, e1, b ) ) / detA;
    if (v < 0) {
        return false;
    }

    intersection.distance = std::numeric_limits<float>::max();

    if ( (u + v) <= 1 ) {
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
