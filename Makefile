CC = g++
CFLAGS= -g -Wall -I. -DGL_GLEXT_PROTOTYPES -Iglut-3.7.6-bin -fopenmp
LDFLAGS= -g -Wall -fopenmp

#Libraries
CCOPTS = -c
LIBRARIES = -lglut -lGL -lGLU

#Final Files and Intermediate .o files
RM = /bin/rm -f
all: main

TARGET = fluid 
main: Point3D.o Particle.o ParticleSystem.o Vector.o SpatialGrid.o window.o
	$(CC) $(LDFLAGS) -o $(TARGET) $^  $(LIBRARIES) 
	$(RM) test

%.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@
	
clean:
	$(RM) test *.o $(TARGET)

rmtest:
	$(RM) test

