all: compile run

compile:
	g++ -o EBGstop -fopenmp EBGstop.cpp

run: 
	./EBGstop
