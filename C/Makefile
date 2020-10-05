all: compile run

compile:
	g++ -o out -I. *.cpp libraylib.a -ldl -lX11 -lpthread -lczmq

run:
	./out