all: paint xrd
xrd: xrd.o domain.o mscs.o rect.o spline.o fft.o
	g++ -o xrd xrd.o domain.o mscs.o rect.o spline.o fft.o
paint: paint.o domain.o mscs.o rect.o spline.o fft.o
	g++ -o paint paint.o domain.o mscs.o rect.o spline.o fft.o
xrd.o: xrd.cpp 
	g++ -c xrd.cpp
paint.o: paint.cpp 
	g++ -c paint.cpp
spline.o: spline.cpp 
	g++ -c spline.cpp
domain.o: domain.cpp
	g++ -c domain.cpp
mscs.o: mscs.cpp
	g++ -c mscs.cpp
rect.o: rect.cpp
	g++ -c rect.cpp
fft.o: fft.cpp
	g++ -c fft.cpp
