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
    int sphereIndex;
};

struct Camera {
    float focalLength;
    vec4  position;
    mat4  rotation;
    float moveSpeed;
    float rotationSpeed;
};

#define SCREEN_WIDTH  1600
#define SCREEN_HEIGHT 1600
#define FULLSCREEN_MODE true
#define SHADOW_BIAS 0.00001f
#define MIN_DEPTH 20
#define RR_PROB 0.70f
#define MAX_SAMPLES 1
#define DROP_FACTOR 1

vector<Triangle> triangles;
vector<Sphere> spheres;
Camera camera;
int max_depth = 0;

vec3 image[SCREEN_WIDTH][SCREEN_HEIGHT];
int iterations = 0;

bool Update();
void Draw(screen* screen);
vec3 castRay(bool insideObject, float &absorbDistance, vec4 &orig, vec4 &dir, int depth);
bool TriangleIntersection(
      vec4 start,
      vec4 dir,
      Triangle triangle,
      Intersection& intersection );
bool SphereIntersection(
      vec4 start,
      vec4 dir,
      Sphere sphere,
      Intersection& intersection );
bool ClosestIntersection(
      vec4 start,
      vec4 dir,
      const vector<Triangle>& triangles,
      const vector<Sphere>& spheres,
      Intersection& closestIntersection );

bool solveQuadratic(const float &a, const float &b, const float &c, float &x0, float &x1);
void RotateX( mat4& rotation, float rad );
void RotateY( mat4& rotation, float rad );
void InitialiseParams(); // Initialise camera and light

// tone mapping methods:
vec3 reinhard(vec3 color);
vec3 filmic(vec3 color);
vec3 hable(vec3 color);

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
        LoadTestModel( triangles, spheres );
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
            float cartita; // Cârtiță
            vec3 pixelColor = castRay(false, cartita, camera.position, primaryRay, 0);
            image[col][row] += pixelColor;
            vec3 tonedColor = filmic(image[col][row] / (float) iterations);
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
    cout << " Max Depth:             " << max_depth << endl;
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
                // case SDLK_w:
                //     // Translate camera forward
                //     camera.position += CS * R * vec4(0, 0, 0.1f, 0);
                //     break;
                // case SDLK_s:
                //     // Translate camera backwards
                //     camera.position -= CS * R * vec4(0, 0, 0.1f, 0);
                //     break;
                // case SDLK_a:
                //     // Translate camera left
                //     camera.position -= CS * R * vec4(0.1f, 0, 0, 0);
                //     break;
                // case SDLK_d:
                //     // Translate camera right
                //     camera.position += CS * R * vec4(0.1f, 0, 0, 0);
                //     break;
                // case SDLK_q:
                //     // Translate camera up
                //     camera.position -= CS * R * vec4(0, 0.1f, 0, 0);
                //     break;
                // case SDLK_e:
                //     // Translate camera down
                //     camera.position += CS * R * vec4(0, 0.1f, 0, 0);
                //     break;
	            // case SDLK_LEFT:
                //     // Rotate camera left
                //     RotateY(camera.rotation, -RS);
		        //     break;
	            // case SDLK_RIGHT:
                //     // Rotate camera right
                //     RotateY(camera.rotation,  RS);
		        //     break;
                // case SDLK_UP:
                //     // Rotate camera up
                //     RotateX(camera.rotation,  RS);
		        //     break;
	            // case SDLK_DOWN:
                //     // Rotate camera down
                //     RotateX(camera.rotation, -RS);
		        //     break;
                // case SDLK_SPACE:
                //     // Reset camera
                //     InitialiseParams();
		        //     break;
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

vec3 castRay(bool insideObject, float &absorbDistance, vec4 &orig, vec4 &dir, int depth)  {
    if (depth > max_depth) {
        max_depth = depth;
    }

    // Russian roulette
    float rr_prob = 1.0f;
    if (depth > MIN_DEPTH) {
        rr_prob = RR_PROB;
        float r_i = distribution(engine);
        if (r_i < 1.0f - rr_prob) {
            return vec3(0,0,0);
        }
    }

    // Compute intersection
    Intersection primaryIntersect;
    if (!ClosestIntersection(orig, dir, triangles, spheres, primaryIntersect)) {
        return vec3(0,0,0);
    }
    // Surface properties
    vec4  normalH;
    vec3  normal;
    int   type;
    float ior;
    vec3  sigma;
    float emission;
    float reflectivity;
    vec3  surfaceColor;

    if (primaryIntersect.triangleIndex > -1) {
        Triangle surface = triangles[primaryIntersect.triangleIndex];
        normalH      = surface.normal;
        normal       = (vec3) normalH;
        type         = surface.isTranparent;
        ior          = surface.IOR;
        sigma        = surface.sigma;
        emission     = surface.emission;
        reflectivity = surface.reflectivity;
        surfaceColor = surface.color;
    } else if (primaryIntersect.sphereIndex > -1) {
        Sphere surface = spheres[primaryIntersect.sphereIndex];
        normal  = (vec3) (primaryIntersect.position - surface.center);
        normal  = normalize(normal);
        normalH = vec4(normal.x, normal.y, normal.z, 1.0f);
        type         = surface.isTranparent;
        ior          = surface.IOR;
        sigma        = surface.sigma;
        emission     = surface.emission;
        reflectivity = surface.reflectivity;
        surfaceColor = surface.color;
    } else {
        cout << "Invalid surface intersection index" << endl;
        exit(1);
    }
    float rnd;

    vec3 shading;
    if (type == TRANSPARENT) {
        /* ----- Transparent BSDF ----- */

        vec3 I = (vec3) dir;
        // n2 -> n1, assume air
        float n1 = 1.000277f;
        float n2 = ior;
        bool enteringObject = true;

        // check if the ray is inside the object
        if (dot(normal, I) > 0) {
            normal  = -normal;
            normalH = -normalH;
            swap(n1, n2);
            enteringObject = false;
        }

        // Schlick's aproximation
        float R0 = (n1 - n2) / (n1 + n2);
              R0 = R0 * R0;
        float cosTheta1 = -dot(normal, I);
        float Rprob = R0 + (1.0f - R0) * pow(1 - cosTheta1, 5);
        float eta = n1/n2;
        float k = 1.0f - eta * eta * (1.0f - cosTheta1 * cosTheta1);
        rnd = distribution(engine);
        if (k >= 0 && rnd > Rprob) {
            // Refraction
            vec3 refractedRay3 = refract(normalize(I), normalize(normal), eta);
            vec4 refractedRay  = vec4(refractedRay3.x,
                                      refractedRay3.y,
                                      refractedRay3.z,
                                      1.0f);
            vec4 refractPoint = primaryIntersect.position - normalH * SHADOW_BIAS;
            shading = castRay(enteringObject, absorbDistance, refractPoint, refractedRay, depth + 1);
            if (enteringObject) {
                // Beer's Law
                vec3 absorbColor = exp((-sigma) * absorbDistance);
                shading *= absorbColor;
                absorbDistance = 0.0f;
            } else {
                // Ray is inside object
                absorbDistance += distance((vec3) orig, (vec3) refractPoint);
            }
        } else {
            // Reflection
            vec3 reflectedRay3 = reflect(normalize(I), normalize(normal));
            vec4 reflectedRay  = vec4(reflectedRay3.x,
                                      reflectedRay3.y,
                                      reflectedRay3.z,
                                      1.0f);
            vec4 reflectPoint = primaryIntersect.position + normalH * SHADOW_BIAS;
            shading = castRay(!enteringObject, absorbDistance, reflectPoint, reflectedRay, depth + 1);
            if (!enteringObject) {
                // Ray is inside object
                absorbDistance += distance((vec3) orig, (vec3) reflectPoint);
            }
        }
        shading /= rr_prob;
    } else if (type == OPAQUE) {
        /* ----- Micro-facet BRDF ----- */

        // Sample from the hemisphere
        vec3 N, Nt, Nb;
        N = normal;
        createCoordinateSystem(N, Nt, Nb);

        vec3 emmittedLight = emission * vec3(1, 1, 1);
        vec4 bouncePoint = primaryIntersect.position + normalH * SHADOW_BIAS;

        float fs = 0.0f;
        if (reflectivity != 0.0f) {
            // Schlick's aproximation
            float n1 = 1.000277f;
            float n2 = ior;
            float R0 = (n1 - n2) / (n1 + n2);
                  R0 = R0 * R0;
            float cosTheta1 = -dot(normal, normalize((vec3) dir));
            float Rprob = R0 + (1.0f - R0) * pow(1 - cosTheta1, 5);
            fs = reflectivity + (1.0f - reflectivity) * Rprob;
        }
        rnd = distribution(engine);

        vec3 indirectLight = vec3(0,0,0);
        if (rnd >= fs) {
            // Spawn a diffuse ray
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

                    indirectLight += t * castRay(insideObject, absorbDistance, bouncePoint, sampleWorld, depth + 1);
                }
            // divide by number of samples, the PDF and the russian roulette factor
            indirectLight /= samples * pdf * rr_prob;
            shading = (indirectLight + emmittedLight) * surfaceColor / (float) M_PI;
        } else {
            // Spawn a specular ray
                vec3 reflectedRay3 = reflect((vec3) dir, normal);
                vec4 reflectedRay  = vec4(reflectedRay3.x,
                                          reflectedRay3.y,
                                          reflectedRay3.z,
                                          1.0f);
                indirectLight = castRay(insideObject, absorbDistance, bouncePoint, reflectedRay, depth + 1);
                indirectLight /= rr_prob;
                shading = indirectLight;
        }

        absorbDistance = (insideObject) ?
                         (absorbDistance + distance((vec3) orig, (vec3) bouncePoint)) :
                          0.0f;
    }

    return shading;
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

bool SphereIntersection( vec4 start, vec4 dir, Sphere sphere, Intersection& intersection ) {
    float t, t0, t1;
    // analytic solution (geometric solution overflows ... )
    vec3  L = (vec3) start - (vec3) sphere.center;
    float a = dot( (vec3) dir, (vec3) dir);
    float b = 2 * dot((vec3) dir, L);
    float c = dot(L, L) - (sphere.radius * sphere.radius);
    if (!solveQuadratic(a, b, c, t0, t1)) return false;

    if (t0 > t1) std::swap(t0, t1);
    if (t0 < 0) {
        t0 = t1; // if t0 is negative, let's use t1 instead
        if (t0 < 0) return false; // both t0 and t1 are negative
    }

    intersection.distance = t0;
    intersection.position = start + dir * t0;
    return true;
}

bool ClosestIntersection(
      vec4 start,
      vec4 dir,
      const vector<Triangle>& triangles,
      const vector<Sphere>& spheres,
      Intersection& closestIntersection ) {

    bool found = false;
    float minDistance = std::numeric_limits<float>::max();
    Intersection localIntersection;

    for( int i = 0; i < triangles.size(); i++ ) {
        if (TriangleIntersection( start, dir, triangles[i], localIntersection )) {
            if (localIntersection.distance < minDistance) {
                found = true;
                closestIntersection = localIntersection;
                closestIntersection.triangleIndex =  i;
                closestIntersection.sphereIndex   = -99;
                minDistance = localIntersection.distance;
            }
        }
    }

    for( int i = 0; i < spheres.size(); i++ ) {
        if (SphereIntersection( start, dir, spheres[i], localIntersection )) {
            if (localIntersection.distance < minDistance) {
                found = true;
                closestIntersection = localIntersection;
                closestIntersection.triangleIndex = -99;
                closestIntersection.sphereIndex   =  i;
                minDistance = localIntersection.distance;
            }
        }
    }

    return found;
}

bool solveQuadratic(const float &a, const float &b, const float &c, float &x0, float &x1) {
    float delta = b * b - 4 * a * c;
    if (delta < 0) {
        return false;
    }
    else if (delta == 0) {
        x0 = x1 = - 0.5 * b / a;
    } else {
        float q = (b > 0) ?
            -0.5 * (b + sqrt(delta)) :
            -0.5 * (b - sqrt(delta));
        x0 = q / a;
        x1 = c / q;
    }
    if (x0 > x1) std::swap(x0, x1);

    return true;
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

vec3 reinhard(vec3 color) {
    vec3 tonedColor;
    tonedColor.x = color.x / (1.0f + color.x);
    tonedColor.y = color.y / (1.0f + color.y);
    tonedColor.z = color.z / (1.0f + color.z);
    return tonedColor;
}

vec3 filmic(vec3 color) {
    float A = 2.51f;
    float B = 0.03f;
    float C = 2.43f;
    float D = 0.59f;
    float E = 0.14f;

    vec3 tonedColor;
    tonedColor.x = ((color.x*(A*color.x+B))/(color.x*(C*color.x+D)+E));
    tonedColor.y = ((color.y*(A*color.y+B))/(color.y*(C*color.y+D)+E));
    tonedColor.z = ((color.z*(A*color.z+B))/(color.z*(C*color.z+D)+E));
    return tonedColor;
}

vec3 hable(vec3 color) {
    float A = 0.15;
    float B = 0.50;
    float C = 0.10;
    float D = 0.20;
    float E = 0.02;
    float F = 0.30;

    vec3 tonedColor;
    tonedColor.x = ((color.x*(A*color.x+C*B)+D*E)/(color.x*(A*color.x+B)+D*F))-E/F;
    tonedColor.y = ((color.y*(A*color.y+C*B)+D*E)/(color.y*(A*color.y+B)+D*F))-E/F;
    tonedColor.z = ((color.z*(A*color.z+C*B)+D*E)/(color.z*(A*color.z+B)+D*F))-E/F;
    tonedColor.x = pow(tonedColor.x, 1.0f / 2.2f);
    tonedColor.y = pow(tonedColor.y, 1.0f / 2.2f);
    tonedColor.z = pow(tonedColor.z, 1.0f / 2.2f);
    return tonedColor;
}
