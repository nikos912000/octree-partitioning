Octree Partitioning 
===================

An implementation for Octree partitioning 3-D points into spatial subvolumes, in C, using PTHREADS and OPENMP.

1st Course Assignment for Parallel and Distributed Computing Systems (2013).

How to use
----------
Just run 'make' command in a unix-based system and run the executable given the appropriate arguments! 

****Arguments for PTHREADS version****


*number of threads, number of 3-D points, maximum number of points inside a cube* 

(eg. ./octree_pthreads 8 100 4)


****Arguments for OPENMP version:**** 

*number of 3-D points, maximum number of points inside a cube* 

(eg. ./octree_openmp 1048576 20)    

Output
------
Info messages in stdout. Results are saved to files.