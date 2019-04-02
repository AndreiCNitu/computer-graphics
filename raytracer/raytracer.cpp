#include <iostream>
#include <glm/glm.hpp>
#include <SDL.h>
#include "SDLauxiliary.h"
#include "TestModelH.h"
#include "ModelLoader.h"
#include <stdint.h>
#include <limits>
#include <math.h>
#include <random>

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

#define SCREEN_WIDTH  1080
#define SCREEN_HEIGHT 1080
#define FULLSCREEN_MODE true
#define SHADOW_BIAS 0.00064f
#define MIN_DEPTH 5
#define RR_PROB 0.95f
#define MAX_SAMPLES 1
#define DROP_FACTOR 1

vector<Triangle> triangles;
Camera camera;

vec3 image[SCREEN_WIDTH][SCREEN_HEIGHT];
int iterations = 0;

bool Update();
void Draw(screen* screen);
vec3 castRay(vec4 &orig, vec4 &dir, int depth);
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
void InitialiseParams(); // Initialise camera and light
void filmic(float &x);
void hable(float &x);

default_random_engine engine(std::random_device{}());
uniform_real_distribution<float> distribution(0.0f, 1.0f);
vec3 uniformSampleHemisphere(const float t, const float p);
void createCoordinateSystem(const vec3 &N, vec3 &Nt, vec3 &Nb);

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
        exit(1);
    }

    SDL_SaveImage( screen, "screenshot.bmp" );

    KillSDL(screen);
    return 0;
}

void Draw(screen* screen) {
    /* Clear buffer */
    memset(screen->buffer, 0, screen->height * screen->width * sizeof(uint32_t));
    iterations++;

    // Loop throught all pixels
    #pragma omp parallel for schedule(dynamic)
    for (int row = 0; row < SCREEN_HEIGHT; row++) {
        for (int col = 0; col < SCREEN_WIDTH; col++) {
            float col_offset = distribution(engine) - 0.5f;
            float row_offset = distribution(engine) - 0.5f;
            vec4 primaryRay = camera.rotation * vec4(col - SCREEN_WIDTH  / 2.0f + col_offset,
                                                     row - SCREEN_HEIGHT / 2.0f + row_offset,
                                                     camera.focalLength, 0.0);

            vec3 pixelColor = castRay(camera.position, primaryRay, 0);
            image[col][row] += pixelColor;
            vec3 tonedColor = image[col][row] / (float) iterations;
            filmic(tonedColor.x);
            filmic(tonedColor.y);
            filmic(tonedColor.z);
            PutPixelSDL(screen, col, row, tonedColor);
        }
    }
}

bool Update() {
    static int t = SDL_GetTicks();
    /* Compute frame time */
    int t2 = SDL_GetTicks();
    float dt = float(t2-t);
    t = t2;
    system("clear");
    cout << "-------------------------------" << endl;
    cout << " Iterations:            " << iterations << endl;
    cout << " Samples / n-th bounce: " << MAX_SAMPLES << " / (" << DROP_FACTOR << " ^ d)" << endl;
    cout << " Frame render time:     " << dt << " ms" << endl;
    cout << "-------------------------------" << endl;

    SDL_Event e;
    while(SDL_PollEvent(&e)) {
        if (e.type == SDL_QUIT) {
	        return false;
	    } else if (e.type == SDL_KEYDOWN) {
            // iterations = 0;
            // for (int y = 0; y < SCREEN_WIDTH; y++)
            //     for (int x = 0; x < SCREEN_HEIGHT; x++)
            //         image[y][x] = vec3(0,0,0);

	        int key_code = e.key.keysym.sym;
            float CS = camera.moveSpeed;
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
                case SDLK_SPACE:
                    // Reset camera
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
}

vec3 castRay(vec4 &orig, vec4 &dir, int depth)  {
    float rr_prob = 1.0f;
    if (depth > MIN_DEPTH) {
        rr_prob = RR_PROB;
        float r_i = distribution(engine);
        if (r_i < 1.0f - rr_prob) {
            return vec3(0,0,0);
        }
    }

    // compute intersection
    Intersection primaryIntersect;
    if (!ClosestIntersection(orig, dir, triangles, primaryIntersect)) {
        return vec3(0,0,0);
    }
    Triangle surface = triangles[primaryIntersect.triangleIndex];

    // Sample from the hemisphere
    vec3 indirectLight = vec3(0,0,0);
    vec3 N, Nt, Nb;
    N = (vec3) surface.normal;
    createCoordinateSystem(N, Nt, Nb);

    vec3 pointColor;
    vec3 emmittedLight = surface.emission * vec3(1, 1, 1);
    if (surface.type == 0) {
        int samples = max( (int) (MAX_SAMPLES / pow(DROP_FACTOR, depth)), 1 );
        float pdf = 1 / (2 * M_PI);
        for (int i = 0; i < samples; ++i) {
            // create sample in world space
            float t = distribution(engine);
            float p = distribution(engine);
            vec3 sample = uniformSampleHemisphere(t, p);
            vec4 sampleWorld = vec4(
                 sample.x * Nb.x + sample.y * N.x + sample.z * Nt.x,
                 sample.x * Nb.y + sample.y * N.y + sample.z * Nt.y,
                 sample.x * Nb.z + sample.y * N.z + sample.z * Nt.z,
                 1.0f);

            vec4 bouncePoint = primaryIntersect.position + sampleWorld * SHADOW_BIAS;
            indirectLight += t * castRay(bouncePoint, sampleWorld, depth + 1);
        }
        // divide by number of samples, the PDF and the russian roulette factor
        indirectLight /= samples * pdf * rr_prob;
        pointColor = (emmittedLight + indirectLight) * surface.color;
    } else if (surface.type == 1) {
        vec3 reflectedRay3 = reflect((vec3) dir, (vec3) surface.normal);
        vec4 reflectedRay  = vec4(reflectedRay3.x,
                                  reflectedRay3.y,
                                  reflectedRay3.z,
                                  1.0f);
        vec4 bouncePoint = primaryIntersect.position + surface.normal * SHADOW_BIAS;
        indirectLight += castRay(bouncePoint, reflectedRay, depth + 1);
        pointColor = emmittedLight + indirectLight;
    }
    return pointColor / (float) M_PI;
}

void createCoordinateSystem(const vec3 &N, vec3 &Nt, vec3 &Nb)  {
    if (std::fabs(N.x) > std::fabs(N.y)) {
        Nt = vec3(N.z, 0, -N.x) / sqrtf(N.x * N.x + N.z * N.z);
    } else {
        Nt = vec3(0, -N.z, N.y) / sqrtf(N.y * N.y + N.z * N.z);
    }
    Nb = cross(N, Nt);
}

vec3 uniformSampleHemisphere(const float t, const float p) {
    // r1 <- [0,1] => 1 - r1 <- [0,1] => we can use cos(theta) = r1
    float phi = 2 * M_PI * p;
    float y = t; // cos(theta)
    // sin(theta) = srt(1 - cos^2(theta))
    float x = sqrtf(1 - t * t) * cosf(phi); // sin(theta) * cos(phi)
    float z = sqrtf(1 - t * t) * sinf(phi); // sin(theta) * sin(phi)
    return vec3(x, y, z);
}

bool TriangleIntersection( vec4 start, vec4 dir, Triangle triangle, Intersection& intersection ) {

    vec3 dir3 = (vec3) dir;
    vec4 v0 = triangle.v0;
    vec4 v1 = triangle.v1;
    vec4 v2 = triangle.v2;

    // v0 + u e1 + v e2 = s + t d
    vec3 e1 = vec3(v1.x - v0.x, v1.y - v0.y, v1.z - v0.z); // e1 = v1 - v0
    vec3 e2 = vec3(v2.x - v0.x, v2.y - v0.y, v2.z - v0.z); // e2 = v2 - v0
    vec3 b  = vec3(start.x - v0.x, start.y - v0.y, start.z - v0.z); // s - v0

    // v0 + ue1 + ve2 = s + td
    mat3 A( -dir3, e1, e2 );
    const float detA  = glm::determinant( A );

    float t = glm::determinant( mat3(   b, e1, e2 ) ) / detA;
    if (t < 0) {
        return false;
    }
    float u = glm::determinant( mat3( -dir3, b, e2 ) ) / detA;
    if (u < 0) {
        return false;
    }
    float v = glm::determinant( mat3( -dir3, e1, b ) ) / detA;
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
    mat4 R = mat4(vec4(1,        0,  0,        0),
                  vec4(0,  cos(rad), sin(rad), 0),
                  vec4(0, -sin(rad), cos(rad), 0),
                  vec4 (0,        0,  0,       1));
    rotation = R * rotation;
}

void RotateY( mat4& rotation, float rad ) {
    mat4 R = mat4(vec4(cos(rad), 0, -sin(rad), 0),
                  vec4(0,        1,  0,        0),
                  vec4(sin(rad), 0,  cos(rad), 0),
                  vec4(0,        0,  0,        1));
    rotation = R * rotation;
}

void hable(float &x) {
    float A = 0.15;
    float B = 0.50;
    float C = 0.10;
    float D = 0.20;
    float E = 0.02;
    float F = 0.30;

    x = ((x*(A*x+C*B)+D*E)/(x*(A*x+B)+D*F))-E/F;
}

void filmic(float &x) {
    float a = 2.51f;
    float b = 0.03f;
    float c = 2.43f;
    float d = 0.59f;
    float e = 0.14f;
    x = ((x*(a*x+b))/(x*(c*x+d)+e));
}
