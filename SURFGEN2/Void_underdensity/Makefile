CC       = gcc -fno-stack-protector
CFLAGS   = -g -O3 -Wall

OBJS=main.o Percolation_pbc_new.o MarchingCube_me_pbc.o cal_genus_imc.o fixup.o funcs.o 

target: main
.phony: clean


main:   $(OBJS)
	$(CC) $(CFLAGS) -o main $(OBJS)  -lm 
	rm -rf *.o

clean:
	rm -rf main *.o *~
