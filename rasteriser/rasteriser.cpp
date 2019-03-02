#include <iostream>
#include <glm/glm.hpp>
#include <SDL.h>
#include "SDLauxiliary.h"
#include "TestModelH.h"
#include <stdint.h>

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
    mat4  rotation;
};

vector<Triangle> triangles;
struct Camera camera;

#define SCREEN_WIDTH  2560
#define SCREEN_HEIGHT 1600
#define FULLSCREEN_MODE true

bool Update();
void Draw(screen* screen);
void VertexShader( const vec4& vertex, ivec2& pixel );
void Interpolate( ivec2 a, ivec2 b, vector<ivec2>& result );
void DrawLineSDL( screen* screen, ivec2 a, ivec2 b, vec3 color );
void DrawPolygonEdges( screen* screen, const vector<vec4>& vertices );

int main( int argc, char* argv[] ) {
    screen *screen = InitializeSDL( SCREEN_WIDTH, SCREEN_HEIGHT, FULLSCREEN_MODE );

    camera.focalLength = SCREEN_HEIGHT;
    camera.position = vec4( 0.0, 0.0, -3.001f, 1.0);
    camera.rotation = mat4(vec4(1.0, 0.0, 0.0, 1.0),
                           vec4(0.0, 1.0, 0.0, 1.0),
                           vec4(0.0, 0.0, 1.0, 1.0),
                           vec4(0.0, 0.0, 0.0, 1.0));

    LoadTestModel( triangles );

    while ( Update()) {
        Draw(screen);
        SDL_Renderframe(screen);
    }

    SDL_SaveImage( screen, "screenshot.bmp" );

    KillSDL(screen);
    return 0;
}

/*Place your drawing here*/
void Draw(screen* screen) {
    /* Clear buffer */
    memset(screen->buffer, 0, screen->height*screen->width*sizeof(uint32_t));

    for(int i = 0; i < triangles.size(); i++) {
        vector<vec4> vertices(3);

        vertices[0] = triangles[i].v0;
        vertices[1] = triangles[i].v1;
        vertices[2] = triangles[i].v2;
        DrawPolygonEdges( screen, vertices );
    }
}

/*Place updates of parameters here*/
bool Update() {
    static int t = SDL_GetTicks();
    /* Compute frame time */
    int t2 = SDL_GetTicks();
    float dt = float(t2-t);
    t = t2;
    cout << "Render time: " << dt << " ms." << endl;

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
            }
        }
    }
    return true;
}

void VertexShader( const vec4& vertex, ivec2& pixel ) {
    float f = camera.focalLength;
    float X = vertex.x;
    float Y = vertex.y;
    float Z = vertex.z - camera.position.z;

    pixel.x = (f * X / Z) + SCREEN_WIDTH / 2;
    pixel.y = (f * Y / Z) + SCREEN_HEIGHT / 2;
}

void Interpolate( ivec2 a, ivec2 b, vector<ivec2>& result ) {
    int N = result.size();
    vec2 step = vec2(b-a) / float(max(N-1,1));
    vec2 current( a );
    for(int i = 0; i < N; i++) {
        result[i] = current;
        current += step;
    }
}

void DrawLineSDL( screen* screen, ivec2 a, ivec2 b, vec3 color ) {
    ivec2 delta = glm::abs( a - b );
    int pixels  = glm::max( delta.x, delta.y ) + 1;
    vector<ivec2> line( pixels );
    Interpolate( a, b, line );
    for(int i = 0; i < pixels; i++) {
        PutPixelSDL(screen, line[i].x, line[i].y, vec3(1,1,1));
    }
}

void DrawPolygonEdges( screen* screen, const vector<vec4>& vertices ) {
    int V = vertices.size();
    // Transform each vertex from 3D world position to 2D image position:
    vector<ivec2> projectedVertices( V );
    for(int i = 0; i < V; i++) {
        VertexShader( vertices[i], projectedVertices[i] );
    }
    // Loop over all vertices and draw the edge from it to the next vertex:
    for(int i = 0; i < V; i++) {
        int j = (i + 1) % V; // The next vertex
        vec3 color( 1, 1, 1 );
        DrawLineSDL( screen, projectedVertices[i], projectedVertices[j], color );
    }
}
