#ifndef MOVEMENT
#define MOVEMENT

#include <glm/glm.hpp>

using glm::vec4;
using glm::mat4;

void ResetParams(vec4& translation, mat4& rotation, vec4& lightPos) {
    translation = vec4( 0.0f, 0.0f, 0.0f, 1.0f);
    rotation = mat4(vec4(1.0, 0.0, 0.0, 0.0),
                           vec4(0.0, 1.0, 0.0, 0.0),
                           vec4(0.0, 0.0, 1.0, 0.0),
                           vec4(0.0, 0.0, 0.0, 1.0));
    lightPos = vec4(0,-0.5,-0.7,1.0);
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

void move( const uint8_t *keyState, int lX, int lY, float CS,
    float LS, float RS, vec4& translation, mat4& rotation, vec4& lightPos) {
    /* Move camera forward */
    if (keyState[SDL_SCANCODE_W]) {
        translation += vec4(0, 0, -CS, 0);
    }
    /* Move camera backwards */
    if (keyState[SDL_SCANCODE_S]) {
        translation += vec4(0, 0, CS, 0);
    }
    /* Move camera left */
    if (keyState[SDL_SCANCODE_A]) {
        translation += vec4(CS, 0, 0, 0);
    }
    /* Move camera right */
    if (keyState[SDL_SCANCODE_D]) {
        translation += vec4(-CS, 0, 0, 0);
    }
    /* Move camera up */
    if (keyState[SDL_SCANCODE_Q]) {
        translation += vec4(0, CS, 0, 0);
    }
    /* Move camera down */
    if (keyState[SDL_SCANCODE_E]) {
        translation += vec4(0, -CS, 0, 0);
    }
    /* Translate light forward */
    if (keyState[SDL_SCANCODE_I]) {
        lightPos += vec4(0, 0, LS, 0);
    }
    /* Translate light backwards */
    if (keyState[SDL_SCANCODE_K]) {
        lightPos += vec4(0, 0, -LS, 0);
    }
    /* Translate light left */
    if (keyState[SDL_SCANCODE_J]) {
        lightPos += vec4(-LS, 0, 0, 0);
    }
    /* Translate light right */
    if (keyState[SDL_SCANCODE_L]) {
        lightPos += vec4(LS, 0, 0, 0);
    }
    /* Translate light up */
    if (keyState[SDL_SCANCODE_U]) {
        lightPos += vec4(0, -LS, 0, 0);
    }
    /* Translate light down */
    if (keyState[SDL_SCANCODE_O]) {
        lightPos += vec4(0, LS, 0, 0);
    }
    /* Reset camera and light */
    if (keyState[SDL_SCANCODE_SPACE]) {
        ResetParams(translation, rotation, lightPos);
    }
    /* Rotate camera left */
    if (lX > 0) {
        RotateY(rotation, -RS);
    }
    /* Rotate camera right */
    if (lX < 0) {
        RotateY(rotation,  RS);
    }
    /* Rotate camera up */
    if (lY > 0) {
        RotateX(rotation,  RS);
    }
    /* Rotate camera down */
    if (lY < 0) {
        RotateX(rotation, -RS);
    }

    /* Rotate camera left */
    if (keyState[SDL_SCANCODE_LEFT]) {
        RotateY(rotation, RS);
    }
    /* Rotate camera right */
    if (keyState[SDL_SCANCODE_RIGHT]) {
        RotateY(rotation, -RS);
    }
    /* Rotate camera up */
    if (keyState[SDL_SCANCODE_UP]) {
        RotateX(rotation, -RS);
    }
    /* Rotate camera down */
    if (keyState[SDL_SCANCODE_DOWN]) {
        RotateX(rotation, RS);
    }

}

#endif
