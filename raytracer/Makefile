FILE=raytracer
GLMDIR=../glm/
EXEC=$(FILE)

# default build settings
CC_OPTS=-c -pipe -Wall -Wno-switch -O3 -qopenmp
LN_OPTS=-qopenmp
CC=icpc

# SDL options
SDL_CFLAGS := $(shell sdl2-config --cflags)
GLM_CFLAGS := -I$(GLMDIR)
SDL_LDFLAGS := $(shell sdl2-config --libs)

# Objects
OBJ = $(FILE).o
$(FILE).o : $(FILE).cpp SDLauxiliary.h TestModelH.h
	$(CC) $(CC_OPTS) -o $(FILE).o $(FILE).cpp $(SDL_CFLAGS) $(GLM_CFLAGS)

# Build rule
build : $(OBJ) Makefile
	$(CC) $(LN_OPTS) -o $(EXEC) $(OBJ) $(SDL_LDFLAGS)

run :
	make clean
	make build
	./raytracer

# Clean rule
.PHONY: clean

clean:
	rm -f $(EXEC) *.o *.bmp
