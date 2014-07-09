CC = gcc
FLG = -O4
NAME = octree_openmp

octree_pthreads: octree_openmp.c octree_openmp.h

	$(CC) octree_openmp.c -fopenmp -lm -o $(NAME)

clean:
	rm -f *.o *.out *.exe
	rm -f *.bin