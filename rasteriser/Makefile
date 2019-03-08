FILE=rasteriser
GLMDIR=../glm/
EXEC=$(FILE)

OS := $(shell uname)

# Build options, OS independent
ifeq ($(OS), Darwin)
CC_OPTS=--std=c++11 -c -pipe -Wall -Wno-switch -Ofast -qopenmp
LN_OPTS=-qopenmp
CC=icpc
else
CC_OPTS=--std=c++11 -c -pipe -Wall -Wno-switch -Ofast -fopenmp
LN_OPTS=-fopenmp
CC=g++
endif

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

run-cornell :
	make clean
	make build
	./$(EXEC) --realtime --cornell-box

run-cube:
	make clean
	make build
	./$(EXEC) --load-box cube.obj

# Clean rule
.PHONY: clean

clean:
	rm -f $(EXEC) *.o *.bmp
