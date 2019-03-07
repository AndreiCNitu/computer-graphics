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
    float rotationSpeed;
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
#define SHADOW_BIAS 0.00064f
#define PHONG_N 10
#define KD 0.60f
#define KS 0.008f
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
vec3 DirectLight(   const Intersection& intersection );
vec3 DiffuseLight(  const Intersection& intersection );
vec3 SpecularLight( const Intersection& intersection );
void InitialiseParams(); // Initialise camera and light

int main( int argc, char* argv[] ) {

    if (argc > 4) {
        cout << "Too many arguments." << endl;
        cout << "Usage: ./raytracer [N renders] --cornell-box or ./raytracer [N renders] --load[-box] <model.obj>" << endl;
        exit(1);
    } else if (argc == 3 && ( strcmp("--cornell-box", argv[2]) == 0 ) ||
               argc == 1) {
        LoadTestModel( triangles );
        cout << "Cornell box test model loaded successfully." << endl;
    } else if (argc == 4 && strcmp("--load",          argv[2]) == 0 ) {
        if (LoadModel( triangles, argv[3] )) {
            cout << argv[3] << " model loaded successfully." << endl;
        } else {
            cout << "Unable to load " << argv[3] << endl;
            exit(1);
        }
    } else if (argc == 4 && strcmp("--load-box",      argv[2]) == 0 ) {
        LoadCornellBox( triangles );
        cout << "Cornell box loaded successfully." << endl;
        vector<Triangle> objTriangles;
        if (LoadModel( objTriangles, argv[3] )) {
            triangles.insert(triangles.end(), objTriangles.begin(), objTriangles.end());
            cout << argv[3] << " model loaded successfully." << endl;
        } else {
            cout << "Unable to load " << argv[3] << endl;
            exit(1);
        }
    } else {
        cout << "Unknown command." << endl;
        cout << "Usage: ./raytracer --cornell-box or ./raytracer --load <model.obj>" << endl;
        exit(1);
    }

    screen *screen = InitializeSDL( SCREEN_WIDTH, SCREEN_HEIGHT, FULLSCREEN_MODE );
    InitialiseParams();

    if ( strcmp("--realtime", argv[1]) == 0 ) {
        while ( Update() ) {
            Draw(screen);
            SDL_Renderframe(screen);
        }
    } else if ( strcmp("--once", argv[1]) == 0 ) {
        Update();
        Draw(screen);
        SDL_Renderframe(screen);
        Update();
    } else {
        cout << "Unknown command." << endl;
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
    #pragma omp parallel for schedule(dynamic)
    for (int row = 0; row < SCREEN_HEIGHT; row++) {
        for (int col = 0; col < SCREEN_WIDTH; col++) {
            vec3 pixelColor = vec3(0,0,0);
            for(int sample = 0; sample < AA_SAMPLES; sample++) {
                vec4 primaryRay = camera.rotation * vec4(col - SCREEN_WIDTH  / 2.0f + col_offset[sample],
                                                         row - SCREEN_HEIGHT / 2.0f + row_offset[sample],
                                                         camera.focalLength, 0.0);
                Intersection primaryIntersect;
                if (ClosestIntersection(camera.position, primaryRay, triangles, primaryIntersect)) {
                    if (triangles[primaryIntersect.triangleIndex].color != vec3(0.15f, 0.15f, 0.75f)) {
                        pixelColor += triangles[primaryIntersect.triangleIndex].color
                                    * (DirectLight(primaryIntersect) + light.indirect);
                    } else {
                        vec3 surfaceColor = triangles[primaryIntersect.triangleIndex].color;
                        pixelColor += DiffuseLight(primaryIntersect ) * KD * surfaceColor
                                    + SpecularLight(primaryIntersect) * KS;
                    }
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
            float CS = camera.moveSpeed;
            float LS = light.moveSpeed;
            float RS = camera.rotationSpeed;
            mat4  R  = camera.rotation;
	        switch(key_code) {
                case SDLK_w:
                    // Translate camera forward
                    camera.position += CS * R * vec4(0, 0, 0.1f, 0);
                    break;
                case SDLK_s:
                    // Translate camera backwards
                    camera.position -= CS * R * vec4(0, 0, 0.1f, 0);
                    break;
                case SDLK_a:
                    // Translate camera left
                    camera.position -= CS * R * vec4(0.1f, 0, 0, 0);
                    break;
                case SDLK_d:
                    // Translate camera right
                    camera.position += CS * R * vec4(0.1f, 0, 0, 0);
                    break;
                case SDLK_q:
                    // Translate camera up
                    camera.position -= CS * R * vec4(0, 0.1f, 0, 0);
                    break;
                case SDLK_e:
                    // Translate camera down
                    camera.position += CS * R * vec4(0, 0.1f, 0, 0);
                    break;
	            case SDLK_LEFT:
                    // Rotate camera left
                    RotateY(camera.rotation, -RS);
		            break;
	            case SDLK_RIGHT:
                    // Rotate camera right
                    RotateY(camera.rotation,  RS);
		            break;
                case SDLK_UP:
                    // Rotate camera up
                    RotateX(camera.rotation,  RS);
		            break;
	            case SDLK_DOWN:
                    // Rotate camera down
                    RotateX(camera.rotation, -RS);
		            break;
                case SDLK_i:
                    // Translate light forward
                    light.position += LS * vec4(0, 0, 1.0f, 0);
                    break;
                case SDLK_k:
                    // Translate light backwards
                    light.position -= LS * vec4(0, 0, 1.0f, 0);
                    break;
                case SDLK_j:
                    // Translate light left
                    light.position -= LS * vec4(1.0f, 0, 0, 0);
                    break;
                case SDLK_l:
                    // Translate light right
                    light.position += LS * vec4(1.0f, 0, 0, 0);
                    break;
                case SDLK_u:
                    // Translate light up
                    light.position -= LS * vec4(0, 1.0f, 0, 0);
                    break;
                case SDLK_o:
                    // Translate light down
                    light.position += LS * vec4(0, 1.0f, 0, 0);
                    break;
                case SDLK_SPACE:
                    // Reset camera and light
                    InitialiseParams();
		            break;
	            case SDLK_ESCAPE:
                    // Quit
		            return false;
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
    camera.moveSpeed     = 1.0f;
    camera.rotationSpeed = 0.1f;
    light.position = vec4( 0.0f, -0.5f, -1.08f, 1.0f );
    light.color = 14.0f * vec3( 1.0, 1.0, 1.0 );
    light.indirect = 0.3f*vec3( 1, 1, 1 );
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

    if ( (u + v) <= 1) {
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

vec3 DiffuseLight( const Intersection& intersection ) {
    return DirectLight(intersection) + light.indirect;
}

vec3 SpecularLight( const Intersection& intersection ) {
/*     R\   |N  /Li
 *       \  |  /
 *        \ | /
 *_________\|/___________
 * R = 2 • (N • L) N - L
*/
    vec3 incident = (vec3) (intersection.position - light.position);
    vec3 normal   = (vec3) (triangles[intersection.triangleIndex].normal);

    vec3 reflected = 2.0f * dot(normal, incident) * normal - incident;
    vec3 viewDir = (vec3) (intersection.position - camera.position);
    reflected = normalize(reflected);
    viewDir   = normalize(viewDir);

    if (DirectLight(intersection) == vec3(0,0,0)) {
        return vec3(0,0,0);
    }
    return light.color * (float) pow(dot(viewDir, reflected), PHONG_N);
}
