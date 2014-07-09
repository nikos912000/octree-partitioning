CC = gcc
FLG = -O4
NAME = octree_pthreads

octree_pthreads: octree_pthreads.c octree_pthreads.h

	$(CC) octree_pthreads.c -lpthread -lm -o $(NAME)

clean:
	rm -f *.o *.out *.exe
	rm -f *.bin