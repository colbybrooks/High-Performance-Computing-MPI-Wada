SOURCES=rk4.cpp rk4.h wadaMPI.cpp pendulum.h pendulum3.cpp #defines files to be compiled
OBJECTS=wadaMPI.o pendulum3.o rk4.o #defines output files
LIBS=-lm  #math library
CXX=mpicxx \-std=c++11 #defines the compiler and c++ version
CXXFLAGS=-O3 -c  #options are optimize and compile

#defines a target for wada 
wada: $(OBJECTS)
	$(CXX) $(OBJECTS) $(LIBS) -o wadaMPI

#dependencies
wadaMPI.o: rk4.h pendulum.h
rk4.o: rk4.h
pendulum3.o: rk4.h pendulum.h 

#cleans wada
clean:
	@rm -rf $(OBJECTS) a.out core wadaMPI