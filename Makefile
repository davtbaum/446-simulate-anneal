build: main.o adj.o util.o
	gcc -g main.o adj.o util.o -o adj
main.o: main.c glob.h
	gcc -g -c main.c
adj.o: adj.c adj.h glob.h
	gcc -g -c adj.c
util.o: util.c
	gcc -g -c util.c
