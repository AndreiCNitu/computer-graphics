#include <iostream>
#include <glm/glm.hpp>
#include <SDL.h>
#include "SDLauxiliary.h"
#include "TestModelH.h"
#include "ModelLoader.h"
#include "Movement.h"
#include <stdint.h>
#include <glm/ext.hpp>

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
    vec2 uv;
};

struct Vertex {
    vec4 position;
    vec2 uv;
};

vector<Triangle> triangles;
Camera camera;
Light light;

#define SCREEN_WIDTH  1600
#define SCREEN_HEIGHT 1600
#define FULLSCREEN_MODE true

#define EDGE_THRESHOLD_MIN 0.0312
#define EDGE_THRESHOLD_MAX 0.125
#define SUBPIXEL_QUALITY 0.75

// Store 1/z for each pixel
float depthBuffer[SCREEN_HEIGHT][SCREEN_WIDTH];
vec3 image[SCREEN_WIDTH][SCREEN_HEIGHT];
bool fxaa = 0;
vec3 fragColor;
float quality[12] = {1, 1, 1, 1, 1, 1.5, 2.0, 2.0, 2.0, 2.0, 4.0, 8.0};

// TEXTURES
SDL_Surface* texture_surf = SDL_LoadBMP("../Models/textures/spot.bmp");
std::vector<glm::vec3>texture_mat;

vec3 pixel_FXAA (int x, int y);
float toLuma(vec3 rgb);
void clean_FXAA(void);
int pixelFXAA ( int i, int j );

void Update();
void Print_Time();
int convert_texture(void);
void FXAA(screen* screen);
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
                                      vec3 reflectance);
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

    cout << "past main" << endl;

    screen *screen = InitializeSDL( SCREEN_WIDTH, SCREEN_HEIGHT, FULLSCREEN_MODE );
    convert_texture();

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

//     Clear Z-buffer
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
        vertices[0].uv = triangles[i].uv0;
        vertices[1].uv = triangles[i].uv1;
        vertices[2].uv = triangles[i].uv2;
        DrawPolygon( screen, vertices, triangleColor, triangleNormal, currentReflectance);
    }
    if (fxaa) {
        FXAA( screen );
    }
}

void FXAA ( screen* screen ) {

    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < SCREEN_HEIGHT; ++i) {
        for (int j = 0; j < SCREEN_WIDTH; ++j) {
            PutPixelSDL( screen, i, j, pixel_FXAA (i, j));
        }
    }
    clean_FXAA();
}

vec3 pixel_FXAA (int i, int j) {

     vec3 colorCenter = image[i][j];
     vec3 fragColor = image[i][j];
     if (0 >= i || 0 >= j || SCREEN_HEIGHT-1 <= i || SCREEN_WIDTH-1 <= j)  return colorCenter;
     vec2 inverseScreenSize(1.0f/SCREEN_WIDTH, 1.0f/SCREEN_HEIGHT);

     float center = toLuma(colorCenter);

     float lumaDown = toLuma(image[i][j-1]);
     float lumaUp = toLuma(image[i][j+1]);
     float lumaLeft = toLuma(image[i-1][j]);
     float lumaRight = toLuma(image[i+1][j]);

     float lumaMin = min(center,min(min(lumaDown,lumaUp),min(lumaLeft,lumaRight)));
     float lumaMax = max(center,max(max(lumaDown,lumaUp),max(lumaLeft,lumaRight)));
     float lumaRange = lumaMax - lumaMin;

     if(lumaRange < max(EDGE_THRESHOLD_MIN,lumaMax*EDGE_THRESHOLD_MAX)){
         fragColor = colorCenter;
         return image[i][j];
     }

     float lumaDownUp = lumaDown + lumaUp;
     float lumaLeftRight = lumaLeft + lumaRight;

     float lumaDownLeft = toLuma(image[i-1][j-1]);
     float lumaUpRight = toLuma(image[i+1][j+1]);
     float lumaUpLeft = toLuma(image[i-1][j+1]);
     float lumaDownRight = toLuma(image[i+1][j-1]);
     float lumaLeftCorners = lumaDownLeft + lumaUpLeft;
     float lumaDownCorners = lumaDownLeft + lumaDownRight;
     float lumaRightCorners = lumaDownRight + lumaUpRight;
     float lumaUpCorners = lumaUpRight + lumaUpLeft;

     bool isHorizontal = (abs(-2.0f*lumaLeft+lumaLeftCorners)+abs(-2.0*center+lumaDownUp )*2.0+abs(-2.0*lumaRight+lumaRightCorners)
         >= abs((-2.0) * lumaUp+lumaUpCorners)+abs(-2.0*center+lumaLeftRight)*2.0+abs(-2.0*lumaDown+lumaDownCorners));

     float luma1 = isHorizontal ? lumaDown : lumaLeft;
     float luma2 = isHorizontal ? lumaUp : lumaRight;
     float gradient1 = luma1 - center;
     float gradient2 = luma2 - center;

     float gradientScaled = 0.25f*max(abs(gradient1),abs(gradient2));
     float stepLength = isHorizontal ? inverseScreenSize.y : inverseScreenSize.x;
     float lumaLocalAverage = 0.0f;

     if(abs(gradient1) < abs(gradient2)) {
         lumaLocalAverage = 0.5*(luma2 + center);
     } else {
         stepLength = - stepLength;
         lumaLocalAverage = 0.5*(luma1 + center);
     }

     vec2 currentUv(i, j);
     if(isHorizontal){
         currentUv.y += stepLength * 0.5;
     } else {
         currentUv.x += stepLength * 0.5;
     }

     vec2 offset = isHorizontal ? vec2(inverseScreenSize.x,0.0) : vec2(0.0,inverseScreenSize.y);
     vec2 uv1 = currentUv - offset;
     vec2 uv2 = currentUv + offset;

     float lumaEnd1 = toLuma(image[(int)uv1.x][(int)uv1.y]);
     float lumaEnd2 = toLuma(image[(int)uv2.x][(int)uv2.y]);
     lumaEnd1 -= lumaLocalAverage;
     lumaEnd2 -= lumaLocalAverage;

     bool reached1 = abs(lumaEnd1) >= gradientScaled;
     bool reached2 = abs(lumaEnd2) >= gradientScaled;
     bool reachedBoth = reached1 && reached2;

     if(!reached1) uv1 -= offset;
     if(!reached2) uv2 += offset;

     if(!reachedBoth) {
         for(int i = 2; i < 12; i++){
             if(!reached1){
                 lumaEnd1 = toLuma(image[(int)uv1.x][(int)uv1.y]);
                 lumaEnd1 = lumaEnd1 - lumaLocalAverage;
             }
             if(!reached2){
                 lumaEnd2 = toLuma(image[(int)uv2.x][(int)uv2.y]);
                 lumaEnd2 = lumaEnd2 - lumaLocalAverage;
             }
             reached1 = abs(lumaEnd1) >= gradientScaled;
             reached2 = abs(lumaEnd2) >= gradientScaled;
             reachedBoth = reached1 && reached2;

             if(!reached1) {
                 uv1 -= offset * quality[i];
             }
             if(!reached2) {
                 uv2 += offset * quality[i];
             }

             if(reachedBoth){
                 break;
             }
         }
     }

     float distance1 = isHorizontal ? (i - uv1.x) : (j - uv1.y);
     float distance2 = isHorizontal ? (uv2.x - i) : (uv2.y - j);
     bool dir = distance1 < distance2;
     float distanceFinal = min(distance1, distance2);
     float edgeThickness = (distance1 + distance2);
     float pixelOffset = - distanceFinal / edgeThickness + 0.5;
     bool iscenterSmaller = center < lumaLocalAverage;
     bool var = ((dir ? lumaEnd1 : lumaEnd2) < 0.0) != iscenterSmaller;
     float finalOffset = var ? pixelOffset : 0.0;
     float avg = (1.0/12.0) * (2.0 * (lumaDownUp + lumaLeftRight) + lumaLeftCorners + lumaRightCorners);

     float tmp = 0.0;
     if (abs(avg - center)/lumaRange < 0.0) {
         tmp = 0.0;
     } else if (abs(avg - center)/lumaRange > 1.0) {
         tmp = 1.0;
     } else {
         tmp = abs(avg - center)/lumaRange;
     }

     float subPixelOffset1 = tmp;
     float subPixelOffset2 = (-2.0 * subPixelOffset1 + 3.0) * subPixelOffset1 * subPixelOffset1;
     float subPixelOffsetFinal = subPixelOffset2 * subPixelOffset2 * SUBPIXEL_QUALITY;

     finalOffset = max(finalOffset,subPixelOffsetFinal);

     vec2 finalUv(i, j);
     if(isHorizontal){
         finalUv.y += finalOffset * stepLength;
     } else {
         finalUv.x += finalOffset * stepLength;
     }

     image[i][j] = (image[(int)finalUv.x][(int)finalUv.y]+ 1.0f * image[i][j]) / 2.0f;
     return image[i][j];
 }

void clean_FXAA (void) {
    vec3 tmp = vec3(0, 0, 0);

    for (int i = 0; i < SCREEN_HEIGHT; ++i) {
        for (int j = 0; j < SCREEN_WIDTH; ++j) {
            image[i][j] = tmp;
        }
    }
}

float toLuma(vec3 rgb) {
    return sqrt(dot(rgb, vec3(0.299, 0.587, 0.114)));
}

int convert_texture(void) {

    int bpp = texture_surf->format->BytesPerPixel;
    /* Here p is the address to the pixel we want to retrieve */

    std::cout << "texture size is: " << texture_surf->w << " and " << texture_surf->h << std::endl;

    std::cout << "size is: " << texture_mat.size() << std::endl;

    for(int i = 0; i < texture_surf->h; ++i) {
        for(int j = 0; j < texture_surf->w; ++j) {

            Uint8 *p = (Uint8 *)texture_surf->pixels + j * texture_surf->pitch + i * bpp;
            Uint32 data = 0;

            switch (bpp) {
                case 1:
                    data = *p;
                    break;
                case 2:
                    data = *(Uint16 *)p;
                    break;
                case 3:
                    if (SDL_BYTEORDER == SDL_BIG_ENDIAN)
                        data = p[0] << 16 | p[1] << 8 | p[2];
                    else
                        data = p[0] | p[1] << 8 | p[2] << 16;
                    break;
                case 4:
                    data = *(Uint32 *)p;
                    break;

                default:
                    return 0;
            }

        SDL_Color rgb;
        SDL_GetRGB(data, texture_surf->format, &rgb.r, &rgb.g, &rgb.b);

        vec3 tmp;
        tmp.x = (float)rgb.r/255.0f;
        tmp.y = (float)rgb.g/255.0f;
        tmp.z = (float)rgb.b/255.0f;
        texture_mat.push_back(tmp);

        }
    }
    return 0;
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

    if (keyState[SDL_SCANCODE_V]) {
       fxaa = !fxaa;
    }
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
        std::cout << std::fixed << dt / step << " ms" << '\t'
                  << (step * 1000) / dt << " fps" << endl;
        counter = 0;
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

    pixel.uv.x = vertex.uv.x;
	pixel.uv.y = vertex.uv.y;
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

        if (fxaa) image[x][y] = color * returnVector;
        if (!fxaa) PutPixelSDL( screen, x, y, color * returnVector);
    }
}

void DrawPolygon( screen* screen, const vector<Vertex>& vertices,
                                        vec3 color,
                                        vec4 normal,
                                        vec3 reflectance) {
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

    vec2 edge0 = vec2(V2.x - V1.x, V2.y - V1.y);
    vec2 edge1 = vec2(V0.x - V2.x, V0.y - V2.y);
    vec2 edge2 = vec2(V1.x - V0.x, V1.y - V0.y);

    for (int row = ymin; row < ymax; row++) {
        for (int col = xmin; col < xmax; col++) {
            Pixel P;
            P.x = col;
            P.y = row;

            float area = edgeFunction(V0, V1, V2);
            float w0   = edgeFunction(V1, V2, P) / area;
            float w1   = edgeFunction(V2, V0, P) / area;
            float w2   = edgeFunction(V0, V1, P) / area;

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
                if (0.0f > V0.uv.x) {
                    PixelShader(screen, P, color, normal, reflectance);
                } else {
                    P.uv = (V0.uv * V0.zinv * w0 +
                            V1.uv * V1.zinv * w1 +
                            V2.uv * V2.zinv * w2)
                            / (P.zinv);

                    int texture_x = P.uv.x * texture_surf->h;
                    int texture_y = P.uv.y * texture_surf->w;

                    vec3 tmp = texture_mat.at(texture_x*texture_surf->w + texture_y);

                    PixelShader(screen, P, tmp, normal, reflectance);
                }
            }
        }
    }
}

float edgeFunction(Pixel &a, Pixel &b, Pixel &p) {
    return ((p.x - a.x) * (b.y - a.y) - (p.y - a.y) * (b.x - a.x));
}
