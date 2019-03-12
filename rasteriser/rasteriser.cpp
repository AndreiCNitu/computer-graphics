#include <iostream>
#include <glm/glm.hpp>
#include <SDL.h>
#include "SDLauxiliary.h"
#include "TestModelH.h"
#include "ModelLoader.h"
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

/* https://www.h-schmidt.net/FloatConverter/IEEE754.html */
#define SCREEN_WIDTH  1440
#define SCREEN_HEIGHT 900

#define CAMERA_MOVEMENT_SPEED 0.015625f
#define CAMERA_ROTATION_SPEED 0.015625f
#define LIGHT_MOVEMENT_SPEED 0.0078125f
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
void Interpolate( Pixel a, Pixel b, vector<Pixel>& result );
void DrawLineSDL( screen* screen, Pixel a, Pixel b, vec3 color,
                                                    vec4 normal,
                                                    vec3 reflectance );
//void DrawPolygonEdges( screen* screen, const vector<Vertex>& vertices );
void ComputePolygonRows( const vector<Pixel>& vertexPixels,
                               vector<Pixel>& leftPixels,
                               vector<Pixel>& rightPixels );
void DrawPolygonRows( screen* screen, vec3 color,
                     const vector<Pixel>& leftPixels,
                     const vector<Pixel>& rightPixels,
                           vec4 normal,
                           vec3 reflectance );
void DrawPolygon( screen* screen, const vector<Vertex>& vertices,
                                        vec3 color,
                                        vec4 normal,
                                        vec3 reflectance );

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
    for( int y=0; y<SCREEN_HEIGHT; ++y )
           for( int x=0; x<SCREEN_WIDTH; ++x )
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

    /* Move camera forward */
    if (keyState[SDL_SCANCODE_W]) {camera.translation += vec4(0, 0, -CAMERA_MOVEMENT_SPEED, 0);}

    /* Move camera backwards */
    if (keyState[SDL_SCANCODE_S]) {camera.translation += vec4(0, 0, CAMERA_MOVEMENT_SPEED, 0);}

    /* Move camera left */
    if (keyState[SDL_SCANCODE_A]) {camera.translation += vec4(CAMERA_MOVEMENT_SPEED, 0, 0, 0);}

    /* Move camera right */
    if (keyState[SDL_SCANCODE_D]) {camera.translation += vec4(-CAMERA_MOVEMENT_SPEED, 0, 0, 0);}

    /* Move camera up */
    if (keyState[SDL_SCANCODE_Q]) {camera.translation += vec4(0, CAMERA_MOVEMENT_SPEED, 0, 0);}

    /* Move camera down */
    if (keyState[SDL_SCANCODE_E]) {camera.translation += vec4(0, -CAMERA_MOVEMENT_SPEED, 0, 0);}

    /* Translate light forward */
    if (keyState[SDL_SCANCODE_I]) {light.position += vec4(0, 0, LIGHT_MOVEMENT_SPEED, 0);}

    /* Translate light backwards */
    if (keyState[SDL_SCANCODE_K]) {light.position += vec4(0, 0, -LIGHT_MOVEMENT_SPEED, 0);}

    /* Translate light left */
    if (keyState[SDL_SCANCODE_J]) {light.position += vec4(-LIGHT_MOVEMENT_SPEED, 0, 0, 0);}

    /* Translate light right */
    if (keyState[SDL_SCANCODE_L]) {light.position += vec4(LIGHT_MOVEMENT_SPEED, 0, 0, 0);}

    /* Translate light up */
    if (keyState[SDL_SCANCODE_U]) {light.position += vec4(0, -LIGHT_MOVEMENT_SPEED, 0, 0);}

    /* Translate light down */
    if (keyState[SDL_SCANCODE_O]) {light.position += vec4(0, LIGHT_MOVEMENT_SPEED, 0, 0);}

    /* Reset camera and light */
    if (keyState[SDL_SCANCODE_SPACE]) {InitialiseParams();}

    /* Rotate camera left */
    if (lastRecordedX < 0) {RotateY(camera.rotation, -CAMERA_ROTATION_SPEED);}

    /* Rotate camera right */
    if (lastRecordedX > 0) {RotateY(camera.rotation,  CAMERA_ROTATION_SPEED);}

    /* Rotate camera up */
    if (lastRecordedY < 0) {RotateX(camera.rotation,  CAMERA_ROTATION_SPEED);}

    /* Rotate camera down */
    if (lastRecordedY > 0) {RotateX(camera.rotation, -CAMERA_ROTATION_SPEED);}
}

void Print_Time() {

    static int t = 0;
    int counter = 0;
    ++counter;

    if (100 == counter) {
    /* Compute frame time */
    int t2 = SDL_GetTicks();
    float dt = float(t2-t);
    t = t2;
    cout << "Render time: " << dt << " ms." << endl;
    counter = 0;
    t = SDL_GetTicks();
    }
}

void InitialiseParams() {
    camera.focalLength = SCREEN_HEIGHT;
    camera.position = vec4( 0.0, 0.0, -3.001f, 1.0);
    camera.translation = vec4( 0.0f, 0.0f, 0.0f, 1.0f);
    camera.rotation = mat4(vec4(1.0, 0.0, 0.0, 0.0),
                           vec4(0.0, 1.0, 0.0, 0.0),
                           vec4(0.0, 0.0, 1.0, 0.0),
                           vec4(0.0, 0.0, 0.0, 1.0));
    camera.moveSpeed = CAMERA_MOVEMENT_SPEED;
    camera.rotationSpeed = CAMERA_ROTATION_SPEED;

    light.position = vec4(0,-0.5,-0.7,1.0);
    light.power = 14.0f * vec3( 1, 1, 1 );
    light.indirectPowerPerArea = 0.5f * vec3( 1, 1, 1 );
    light.moveSpeed = LIGHT_MOVEMENT_SPEED;
}

void RotateX( mat4& rotation, float rad ) {
    vec4 c1(1,        0,  0,        0);
    vec4 c2(0,  cos(rad), sin(rad), 0);
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

void Interpolate( Pixel a, Pixel b, vector<Pixel>& result ) {

    int N = result.size();
    vec2  step = vec2(b.x - a.x, b.y - a.y) / float(max(N - 1, 1));
    vec2  current = vec2( a.x, a.y );
    for(int i = 0; i < N; i++) {
        result[i].x = round(current.x);
        result[i].y = round(current.y);
        current.x += step.x;
        current.y += step.y;
    }

    float depthStep = (b.zinv - a.zinv) / float(max(N - 1, 1));
	float currentDepth = a.zinv;
    for(int i = 0; i < N; i++) {
        result[i].zinv = currentDepth;
        currentDepth += depthStep;
    }

    vec4 posStep = (b.pos3d * b.zinv - a.pos3d * a.zinv) / float(max(N - 1, 1));
	vec4 currentPos = a.pos3d * a.zinv;
	for( int i = 0; i < N; i++ ) {
		result[i].pos3d = currentPos / result[i].zinv;
		currentPos += posStep;
	}
}

void DrawLineSDL( screen* screen, Pixel a, Pixel b, vec3 color,
                                                    vec4 normal,
                                                    vec3 reflectance ) {
    ivec2 delta = glm::abs( vec2(a.x - b.x, a.y - b.y) );
    int  pixels = glm::max( delta.x, delta.y ) + 1;
    vector<Pixel> line( pixels );
    Interpolate( a, b, line );

    for(int i = 0; i < pixels; i++) {
        if (line[i].x >= 0            &&
            line[i].y >= 0            &&
            line[i].x <  SCREEN_WIDTH &&
            line[i].y <  SCREEN_HEIGHT) {
                PixelShader(screen, line[i], color, normal, reflectance);
        }
    }
}

// void DrawPolygonEdges( screen* screen, const vector<Vertex>& vertices ) {
//     int V = vertices.size();
//     // Transform each vertex from 3D world position to 2D image position:
//     vector<Pixel> projectedVertices( V );
//     for(int i = 0; i < V; i++) {
//         VertexShader( vertices[i], projectedVertices[i] );
//     }
//     // Loop over all vertices and draw the edge from it to the next vertex:
//     for(int i = 0; i < V; i++) {
//         int j = (i + 1) % V; // The next vertex
//         vec3 color( 1, 1, 1 );
//         DrawLineSDL( screen, projectedVertices[i], projectedVertices[j], color );
//     }
// }

void ComputePolygonRows( const vector<Pixel>& vertexPixels,
                               vector<Pixel>& leftPixels,
                               vector<Pixel>& rightPixels ) {
    int maxY = max( vertexPixels[0].y, max(vertexPixels[1].y, vertexPixels[2].y));
    int minY = min( vertexPixels[0].y, min(vertexPixels[1].y, vertexPixels[2].y));
    int rows = maxY - minY + 1;

    leftPixels.resize(rows);
	rightPixels.resize(rows);
    for( int i = 0; i < rows; i++ ) {
		leftPixels[i].x  = +numeric_limits<int>::max();
		rightPixels[i].x = -numeric_limits<int>::max();
        // Initialise y values
        leftPixels[i].y  = i + minY;
        rightPixels[i].y = i + minY;
	}

    int E = vertexPixels.size();
    for( int e = 0; e < E; e++ ) {
        Pixel a = vertexPixels[e];
        Pixel b = vertexPixels[(e + 1) % E];
        ivec2 delta = glm::abs( vec2(a.x - b.x, a.y - b.y) );
        int  pixels = glm::max( delta.x, delta.y ) + 1;
        vector<Pixel> edge(pixels);
        Interpolate( a, b, edge );
        for(int p = 0; p < edge.size(); p++) {
            // relative index
            int idx = edge[p].y - minY;
            if( edge[p].x < leftPixels[idx].x ) {
                leftPixels[idx].x = edge[p].x;
                leftPixels[idx].zinv = edge[p].zinv;
                leftPixels[idx].pos3d = edge[p].pos3d;
            }
            if( edge[p].x > rightPixels[idx].x ) {
                rightPixels[idx].x = edge[p].x;
                rightPixels[idx].zinv = edge[p].zinv;
                rightPixels[idx].pos3d = edge[p].pos3d;
            }
        }
    }
}

void DrawPolygonRows( screen* screen, vec3 color,
                      const vector<Pixel>& leftPixels,
                      const vector<Pixel>& rightPixels,
                            vec4 normal,
                            vec3 reflectance ) {
    for( int row = 0; row < leftPixels.size(); row++ ) {
        DrawLineSDL( screen, leftPixels[row], rightPixels[row], color, normal, reflectance );
    }
}

void DrawPolygon( screen* screen, const vector<Vertex>& vertices,
                                        vec3 color,
                                        vec4 normal,
                                        vec3 reflectance ) {
    int V = vertices.size();
    vector<Pixel> vertexPixels( V );
    for( int i = 0; i < V; i++ ) {
        VertexShader( vertices[i], vertexPixels[i] );
    }
    vector<Pixel> leftPixels;
    vector<Pixel> rightPixels;
    ComputePolygonRows( vertexPixels, leftPixels, rightPixels );
    DrawPolygonRows( screen, color, leftPixels, rightPixels, normal, reflectance );
}
