FILE=raytracer
GLMDIR=../glm/

########
#   Output
EXEC=$(FILE)

# default build settings
CC_OPTS=-c -pipe -Wall -Wno-switch -ggdb -g3 
LN_OPTS=
CC=g++

########
#       SDL options
SDL_CFLAGS := $(shell sdl2-config --cflags)
GLM_CFLAGS := -I$(GLMDIR)
SDL_LDFLAGS := $(shell sdl2-config --libs)

########
#   This is the default action
all:Build


########
#   Object list
#
OBJ = $(FILE).o


########
#   Objects
$(FILE).o : $(FILE).cpp SDLauxiliary.h TestModelH.h
	$(CC) $(CC_OPTS) -o $(FILE).o $(FILE).cpp $(SDL_CFLAGS) $(GLM_CFLAGS)


########
#   Main build rule     
Build : $(OBJ) Makefile
	$(CC) $(LN_OPTS) -o $(EXEC) $(OBJ) $(SDL_LDFLAGS)

clean:
	rm raytracer
	rm raytracer.o
