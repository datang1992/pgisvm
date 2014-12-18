make:	*.o
*.o:	*.cpp
	g++-mp-4.8 -Wall -I/opt/local/include -c *.cpp -fopenmp
	g++-mp-4.8 -L/opt/local/lib *.o -lgsl -lgslcblas -lm -fopenmp -o main
