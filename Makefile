SOURCES=rk4.cpp rk4.h mpi.cpp pendulum.h pendulum.cpp #defines files to be compiled
OBJECTS=mpi.o pendulum.o rk4.o #defines output files
LIBS=-lm  #math library
CXX=mpicxx \-std=c++11 #defines the compiler and c++ version
CXXFLAGS=-O3 -c  #options are optimize and compile

#defines a target for wada 
wada: $(OBJECTS)
	$(CXX) $(OBJECTS) $(LIBS) -o mpi

#dependencies
mpi.o: rk4.h pendulum.h
rk4.o: rk4.h
pendulum.o: rk4.h pendulum.h 

#cleans wada
clean:
	@rm -rf $(OBJECTS) a.out core mpi