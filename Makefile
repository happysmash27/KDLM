all: main.c
	gcc main.c -lm -lgmp -Wall -Wextra -pedantic -o main
