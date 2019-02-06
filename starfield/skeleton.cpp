#include <iostream>
#include <glm/glm.hpp>
#include <SDL.h>
#include "SDLauxiliary.h"
#include "TestModel.h"
#include <stdint.h>

using namespace std;
using glm::vec3;
using glm::mat3;

#define SCREEN_WIDTH 1280
#define SCREEN_HEIGHT 800
#define FULLSCREEN_MODE true

int t;
vector<vec3> leftSide( SCREEN_HEIGHT );
vector<vec3> rightSide( SCREEN_HEIGHT );
vector<vec3> stars( 1000 );

void Update();
void Draw( screen* screen );
void Interpolate( float a, float b, vector<float>& result );
void Interpolate( vec3 a, vec3 b, vector<vec3>& result );

int main( int argc, char* argv[] ) {

    vector<vec3> result(SCREEN_HEIGHT);

    // vec3 topLeft(1,0,0); // red
    // vec3 topRight(0,0,1); // blue
    // vec3 bottomLeft(1,1,0); // yellow
    // vec3 bottomRight(0,1,0); // green

    //Interpolate( topLeft, bottomLeft, leftSide );
    //Interpolate( topRight, bottomRight, rightSide );

    for(int s = 0; s < 1000; s++) {
        stars[s].x = (2*float(rand()) / float(RAND_MAX)) - 1;
        stars[s].y = (2*float(rand()) / float(RAND_MAX)) - 1;
        stars[s].z = float(rand()) / float(RAND_MAX);
    }

    screen *screen = InitializeSDL(SCREEN_WIDTH, SCREEN_HEIGHT, FULLSCREEN_MODE);
    //Set start value for timer.
    t = SDL_GetTicks();

    while(NoQuitMessageSDL()) {
        Draw(screen);
        Update();
        SDL_Renderframe(screen);
    }

    SDL_SaveImage(screen, "screenshot.bmp");
    KillSDL(screen);
    return 0;
}

// Place your drawing here
void Draw(screen* screen) {
    // Clear buffer
    memset(screen->buffer, 0, screen->height*screen->width*sizeof(uint32_t));

    int f = SCREEN_HEIGHT / 2;
    //vec3 colour(1.0, 1.0, 1.0);
    for( int s = 0; s < stars.size(); s++ ) {
        vec3 colour = 0.2f * vec3(1,1,1) / (stars[s].z * stars[s].z);
        int u = f * stars[s].x / stars[s].z + SCREEN_WIDTH/2;
        int v = f * stars[s].y / stars[s].z + SCREEN_HEIGHT/2;
        PutPixelSDL(screen, u, v, colour);
    }

    // for (int row = 0; row < SCREEN_HEIGHT; row++) {
    //     vector<vec3> one_row(SCREEN_WIDTH);
    //     Interpolate(leftSide[row], rightSide[row], one_row);
    //     for (int col = 0; col < SCREEN_WIDTH; col++) {
    //         PutPixelSDL(screen, col, row, one_row[col]);
    //     }
    // }
}

// Place updates of parameters here
void Update() {
    static int t  = SDL_GetTicks();

    // Compute frame time
    int t2 = SDL_GetTicks();
    float dt = float( t2 - t );
    t = t2;

    static float v = 0.0005f;

    for( int s = 0; s < stars.size(); s++ ) {
        stars[s].z -= v * dt;
    }

    for( int s=0; s<stars.size(); s++ ) {
       // Add code for update of stars
       if( stars[s].z <= 0 ) {
           stars[s].z += 1;
       } else if( stars[s].z > 1 ) {
           stars[s].z -= 1;
       }
   }
}

// Interpolate 1
void Interpolate( float a, float b, vector<float>& result ) {
    if( result.size() == 1 ) {
        result[0] = (a + b) / 2;
        return;
    }
    for( int i = 0; i < result.size(); i++ ) {
        result[i] = a + ( ((b - a) * i) / (result.size()-1) );
    }
}

//Interpolate 2
void Interpolate( vec3 a, vec3 b, vector<vec3>& result ) {
    vector<float> xs( result.size() );
    vector<float> ys( result.size() );
    vector<float> zs( result.size() );

    Interpolate( a.x, b.x, xs );
    Interpolate( a.y, b.y, ys );
    Interpolate( a.z, b.z, zs );

    for (int i = 0; i < result.size(); i++) {
        result[i].x = xs[i];
        result[i].y = ys[i];
        result[i].z = zs[i];
    }
}
