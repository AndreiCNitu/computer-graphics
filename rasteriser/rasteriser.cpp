#include <iostream>
#include <glm/glm.hpp>
#include <SDL.h>
#include "SDLauxiliary.h"
#include "TestModelH.h"
#include "ModelLoader.h"
#include "Movement.h"
#include <stdint.h>
#include <omp.h>

using namespace std;
using glm::ivec2;
using glm::vec2;
using glm::vec3;
using glm::mat3;
using glm::vec4;
using glm::mat4;

struct Camera {
    float focalLength;
    vec4  position;
    vec4  translation;
    mat4  rotation;
    float moveSpeed;
    float rotationSpeed;
};

struct Light {
    vec4 position;
    vec3 power;
    vec3 indirectPowerPerArea;
    float moveSpeed;
};

struct Pixel {
    int x;
    int y;
    float zinv;
    vec4 pos3d;
};

struct Vertex {
    vec4 position;
};

vector<Triangle> triangles;
Camera camera;
Light light;

#define SCREEN_WIDTH  1280
#define SCREEN_HEIGHT 800
#define FULLSCREEN_MODE true

// Store 1/z for each pixel
float depthBuffer[SCREEN_HEIGHT][SCREEN_WIDTH];

void Update();
void Print_Time();
void Draw(screen* screen);
void InitialiseParams(); // Initialise camera and light
void RotateX( mat4& rotation, float rad );
void RotateY( mat4& rotation, float rad );
void VertexShader( const Vertex& vertex, Pixel& pixel );
void PixelShader( screen* screen, const Pixel& pixel, vec3 color,
                                                      vec4 normalH,
                                                      vec3 reflectance );
void DrawPolygon( screen* screen, const vector<Vertex>& vertices,
                                        vec3 color,
                                        vec4 normal,
                                        vec3 reflectance );
float edgeFunction(Pixel &a, Pixel &b, Pixel &p);

int main( int argc, char* argv[] ) {
    if (argc > 4) {
        cout << "Too many arguments." << endl;
        cout << "Usage: ./rasteriser [N renders] --cornell-box or ./rasteriser [N renders] --load[-box] <model.obj>" << endl;
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
        cout << "Usage: ./rasteriser --cornell-box or ./rasteriser --load <model.obj>" << endl;
        exit(1);
    }

    screen *screen = InitializeSDL( SCREEN_WIDTH, SCREEN_HEIGHT, FULLSCREEN_MODE );
    InitialiseParams();

    if ( strcmp("--realtime", argv[1]) == 0 ) {
        while (NoQuitMessageSDL()){
            Update();
            Draw(screen);
            SDL_Renderframe(screen);
            Print_Time();
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

/*Place your drawing here*/
void Draw(screen* screen) {
    /* Clear buffer */
    memset(screen->buffer, 0, screen->height * screen->width * sizeof(uint32_t));

    // Clear Z-buffer
    for( int y = 0; y < SCREEN_HEIGHT; ++y )
           for( int x = 0; x < SCREEN_WIDTH; ++x )
            depthBuffer[y][x] = 0.0f;

    #pragma omp parallel for schedule(dynamic)
    for(int i = 0; i < triangles.size(); i++) {
        vec3 triangleColor  = triangles[i].color;
        vec4 triangleNormal = triangles[i].normal;
        vec3 currentReflectance = vec3(1, 1, 1);
        vector<Vertex> vertices(3);

        vertices[0].position = triangles[i].v0;
        vertices[1].position = triangles[i].v1;
        vertices[2].position = triangles[i].v2;
        DrawPolygon( screen, vertices, triangleColor, triangleNormal, currentReflectance );
    }
}

void Update() {
    const uint8_t *keyState = SDL_GetKeyboardState(NULL);
    int lastRecordedX = 0;
    int lastRecordedY = 0;
    SDL_GetRelativeMouseState(&lastRecordedX, &lastRecordedY);

    float CS = camera.moveSpeed;
    float LS = light.moveSpeed;
    float RS = camera.rotationSpeed;

    move(keyState,
         lastRecordedX,
         lastRecordedY,
         CS, LS, RS,
         camera.translation,
         camera.rotation,
         light.position);
}

void Print_Time() {
    static int t = 0;
    static int counter = 0;
    ++counter;
    int step = 100;

    if (step == counter) {
        /* Compute frame time */
        int t2 = SDL_GetTicks();
        float dt = float(t2-t);
        t = t2;
        std::cout.precision(2);
        std::cout << std::fixed << dt / step << " ms" << '\t' << (step * 1000) / dt << " fps" << endl;
        counter = 0;
    }
}

void InitialiseParams() {
    /* https://www.h-schmidt.net/FloatConverter/IEEE754.html */
    camera.focalLength = SCREEN_HEIGHT;
    camera.position = vec4( 0.0, 0.0, -3.001f, 1.0);
    camera.translation = vec4( 0.0f, 0.0f, 0.0f, 1.0f);
    camera.rotation = mat4(vec4(1.0, 0.0, 0.0, 0.0),
                           vec4(0.0, 1.0, 0.0, 0.0),
                           vec4(0.0, 0.0, 1.0, 0.0),
                           vec4(0.0, 0.0, 0.0, 1.0));
    camera.moveSpeed = 0.015625f;
    camera.rotationSpeed = 0.015625f;

    light.position = vec4(0,-0.5,-0.7,1.0);
    light.power = 14.0f * vec3( 1, 1, 1 );
    light.indirectPowerPerArea = 0.5f * vec3( 1, 1, 1 );
    light.moveSpeed = 0.0078125f;
}

mat4 TransformationMatrix() {
    mat4 translationMatrix = mat4(vec4(1.0, 0.0, 0.0, 0.0),
                                  vec4(0.0, 1.0, 0.0, 0.0),
                                  vec4(0.0, 0.0, 1.0, 0.0),
                                  camera.translation);

    mat4 rotationMatrix = camera.rotation;

    mat4 M = mat4(vec4(1.0, 0.0, 0.0, 0.0),
                  vec4(0.0, 1.0, 0.0, 0.0),
                  vec4(0.0, 0.0, 1.0, 0.0),
                  vec4(-camera.position.x, -camera.position.y, -camera.position.z, 1.0f ));

    return translationMatrix * rotationMatrix * M;
}

void VertexShader( const Vertex& vertex, Pixel& pixel ) {

    vec4 projection = TransformationMatrix() * vertex.position;
    float f = camera.focalLength;
    float X = projection.x;
    float Y = projection.y;
    float Z = projection.z;
    if (Z == 0) return;

    pixel.x = (f * X / Z) + SCREEN_WIDTH  / 2;
    pixel.y = (f * Y / Z) + SCREEN_HEIGHT / 2;
    pixel.zinv = 1.0f / Z;

    pixel.pos3d.x = vertex.position.x;
	pixel.pos3d.y = vertex.position.y;
	pixel.pos3d.z = vertex.position.z;
}

void PixelShader( screen* screen, const Pixel& pixel, vec3 color,
                                                      vec4 normalH,
                                                      vec3 reflectance ) {
    int x = pixel.x;
    int y = pixel.y;
    if (x < 0             ||
        y < 0             ||
        x >= SCREEN_WIDTH ||
        y >= SCREEN_HEIGHT) {
            return;
    }

    if( pixel.zinv > depthBuffer[y][x] ) {
        depthBuffer[y][x] = pixel.zinv;

        vec3 lightVector = (vec3) (light.position - pixel.pos3d);
        float radius = length( lightVector );
        lightVector = normalize( lightVector );
        vec4 lightVectorH = vec4(lightVector.x, lightVector.y, lightVector.z, 1.0f);
        vec3 normal = (vec3) normalH;
        vec3 returnVector = light.power * max(dot( lightVector, normal ), 0.0f)
                          / ( (float) (4 * M_PI) * radius * radius );
        returnVector = reflectance * (returnVector + light.indirectPowerPerArea);

        PutPixelSDL( screen, x, y, color * returnVector);
    }
}

void DrawPolygon( screen* screen, const vector<Vertex>& vertices,
                                        vec3 color,
                                        vec4 normal,
                                        vec3 reflectance ) {
    int V = vertices.size();
    vector<Pixel> vertexPixels( V );
    int xmin = numeric_limits<int>::max();
    int ymin = numeric_limits<int>::max();
    int xmax = numeric_limits<int>::min();
    int ymax = numeric_limits<int>::min();

    for( int i = 0; i < V; i++ ) {
        VertexShader( vertices[i], vertexPixels[i] );
        if (vertexPixels[i].x < xmin) {
            xmin = vertexPixels[i].x;
        }
        if (vertexPixels[i].x > xmax) {
            xmax = vertexPixels[i].x;
        }
        if (vertexPixels[i].y < ymin) {
            ymin = vertexPixels[i].y;
        }
        if (vertexPixels[i].y > ymax) {
            ymax = vertexPixels[i].y;
        }
    }

    Pixel V0 = vertexPixels[0];
    Pixel V1 = vertexPixels[1];
    Pixel V2 = vertexPixels[2];
    for (int row = ymin; row < ymax; row++) {
        for (int col = xmin; col < xmax; col++) {
            Pixel P;
            P.x = col;
            P.y = row;

            float area = edgeFunction(V0, V1, V2);
            float w0   = edgeFunction(V1, V2, P) / area;
            float w1   = edgeFunction(V2, V0, P) / area;
            float w2   = edgeFunction(V0, V1, P) / area;

            vec2 edge0 = vec2(V2.x - V1.x, V2.y - V1.y);
            vec2 edge1 = vec2(V0.x - V2.x, V0.y - V2.y);
            vec2 edge2 = vec2(V1.x - V0.x, V1.y - V0.y);

            bool overlaps = true;
            // If the point is on the edge, test if it is a top or left edge,
            // otherwise test if the edge function is positive
            overlaps &= (w0 == 0 ? ((edge0.y == 0 && edge0.x > 0) || (edge0.y < 0)) : (w0 > 0));
            overlaps &= (w1 == 0 ? ((edge1.y == 0 && edge1.x > 0) || (edge1.y < 0)) : (w1 > 0));
            overlaps &= (w2 == 0 ? ((edge2.y == 0 && edge2.x > 0) || (edge2.y < 0)) : (w2 > 0));

            if (overlaps) {
                P.zinv = V0.zinv * w0 + V1.zinv * w1 + V2.zinv * w2;
                P.pos3d = (V0.pos3d * V0.zinv * w0 +
                           V1.pos3d * V1.zinv * w1 +
                           V2.pos3d * V2.zinv * w2)
                         / (P.zinv);
                PixelShader(screen, P, color, normal, reflectance);
            }
        }
    }
}

float edgeFunction(Pixel &a, Pixel &b, Pixel &p) {
    return ((p.x - a.x) * (b.y - a.y) - (p.y - a.y) * (b.x - a.x));
}
