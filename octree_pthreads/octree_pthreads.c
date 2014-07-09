/********************************************************************/
/*    Octree partitioning 3-D points into spatial subvolumes        */
/*    using PTHREADS project 2013                                   */
/*                                                                  */
/*    Implemented by Nikos Katirtzis (nikos912000)                  */
/********************************************************************/


/******************** Includes - Defines ****************/
#include "octree_pthreads.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <pthread.h>
#include <sys/time.h>
#include <time.h>


/*************** Defines - Initializations **************/
//Maximum number of running threads
int num_threads;
//Total number of points
int N;
//Maximum number of points allowed in a (sub)cube
int S;

// Define a matrix for boxes and leaves
Box *box;
Box *leaf; 
double **A, **B;

// Counters
int box_counter = 1; 
int leaf_counter = 0; 
int num_points = 0;
int running_threads = 0;

// Maximum number of levels - to be printed
int num_levels = 0;
// Number of boxes in each level
int *level_boxes;

// Mutexes
pthread_mutex_t lock_leaf = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t lock_box = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t lock_thread = PTHREAD_MUTEX_INITIALIZER;

/******************** Functions ********************/
void* cubeCheck(void *arg);
void* cubeDivision(void *arg);
void searchNeighbours();
void checkB();
void checkBoundaries();

// Write to files
void AFile();
void BFile();
void neighboursFile();



/******************** Main function ********************/
int main(int argc, char** argv)
{
	if(argc < 4)
	{
		printf("Error in arguments! Three arguments required: num_threads, N and S\n");
		return 0;
	}

    int i, j, position, total_time,B_counter = 0;
    
	// get arguments 
	num_threads = atoi(argv[1]);
	N = atoi(argv[2]);
	S = atoi(argv[3]);
	
    /*****************************************************************/
    /* The original N*3 array of random points that belong to the    */
    /* first octant of the unit sphere (x^2+y^=^2+z^2=1, 0<=x,y,z<=1)*/
    /*****************************************************************/
	
    A = (double**)malloc(N * sizeof(double*));
    if (A == NULL)
    {
        exit(1);
    }
    
    // number generator seed for different random sequence per run
    srand(time(NULL));
    // or this command if we want the same sequence for every iteration
    //srand(0);
    
    // all values are double because for large N floats aren't large enough
    double y_max = 0;
    
    for (i = 0; i < N; i++)
    {
        A[i] = (double*)malloc(3 * sizeof(double));
        // points for which it states x^2+y^2+z^2=1 and 0<=x,y,z<=1
        A[i][0] = ((double)rand() / (double)RAND_MAX);
        // y_max = 1 - x^2
        y_max = sqrt(1 - pow(A[i][0], 2));
        A[i][1] = ((double)rand() / (double)RAND_MAX) * y_max;
        // z = 1 - x^2 - y^2
        A[i][2] = sqrt(1 - pow(A[i][0], 2) - pow(A[i][1], 2));  
    }
    
    printf("Generation of points completed!\n");
    
    // start timer
    gettimeofday(&startwtime, NULL);
    
    // memory allocation for array B
    B = (double**) malloc(sizeof(double*) * N);
    if (B == NULL)
    {
        exit(1);
    } 
    
    for (i = 0; i < N; i++)
    {
        B[i] = (double*)malloc(sizeof(double) * 3);
    }
     
    // memory allocation (1 element) for array box
    box = (Box*)malloc(1 * sizeof(Box)); 
    if (box == NULL)
    {
        exit(1);
    }
    
    // initialize 1st box (unit cube))
    box[0].level = 0;
    box[0].boxid = 1;
    box[0].parent = 0;
    box[0].length = 1;
    box[0].center[0] = 0.5;
    box[0].center[1] = 0.5;
    box[0].center[2] = 0.5;
    box[0].start = 0;
    box[0].n = N; 
    box[0].points = (int*)malloc(N * sizeof(int));
    
    for (i = 0; i < 26; i++)
    {
        box[0].colleague[i] = 0;
    }
    
    for (i = 0; i < N; i++)
    {   
        box[0].points[i] = i;
    }
    
    // begin calculations
    position = 0;
    cubeCheck(&position);
    printf("Creation of octree completed!\n");
    
    // find number of cubes in each level
    level_boxes = (int*)malloc((num_levels + 1) * sizeof(int));
    if (level_boxes == NULL){
        exit(1);
    }
    
    printf("Maximum number of levels = %d\n", num_levels);
    printf("Total number of cubes = %d\n", box_counter);
    
    searchNeighbours();
    printf("All colleagues have been found!\n");
    
    //Copy all points from leafs to array B
    for ( i = 0; i < leaf_counter; i++)
    {
        leaf[i].start = B_counter;
        for (j = 0; j < leaf[i].n; j++)
        {
            B[B_counter][0] = A[leaf[i].points[j]][0];
            B[B_counter][1] = A[leaf[i].points[j]][1];
            B[B_counter][2] = A[leaf[i].points[j]][2];
            B_counter++;
        }
    }
    printf("Array B updated!\n");
    /*checkB();
    checkBoundaries();*/
    
    // file insertion
    /*AFile();
    BFile();
    neighboursFile();
    printf("File insertion completed!\n");*/
       
    // stop timer
    gettimeofday(&endwtime, NULL);
    total_time = ((endwtime.tv_sec * 1000000 + endwtime.tv_usec) -(startwtime.tv_sec * 1000000 + startwtime.tv_usec));
    
    printf("Total calculation time is: %d us\n", total_time);
    printf("\nTask completed!\n");
    return (EXIT_SUCCESS);
}


void* cubeCheck(void *arg)
{
    int i;
    int boxIndex = *(int *)arg; 
    Box temp_box, temp_parent;
    
    pthread_mutex_lock(&lock_box);
    temp_box = box[boxIndex];
    temp_parent = box[temp_box.parent - 1];
    pthread_mutex_unlock(&lock_box);

    /*Array with points indexes that belong to the cube*/
    if(temp_box.boxid != 1)
    {
        temp_box.points = (int*)malloc(temp_parent.n * sizeof(int));
        /*Checking how many points are included in the cube*/
        for (i = 0; i < temp_parent.n; i++)
        {   
            if (fabs(temp_box.center[0] - A[temp_parent.points[i]][0]) < temp_box.length / 2)
            {
                if (fabs(temp_box.center[1] - A[temp_parent.points[i]][1]) < temp_box.length / 2)
                {
                    if (fabs(temp_box.center[2] - A[temp_parent.points[i]][2]) < temp_box.length / 2)
                    {      
                        temp_box.n++;
                        temp_box.points[temp_box.n-1] = temp_parent.points[i];
                    }
                }
            }
        }
    } 
    
    if (temp_box.n == 0)
    { 
        // cube has no points (empty)...set boxid = 0 and this child of parent = 0
        pthread_mutex_lock(&lock_box);
        temp_box.boxid = 0;
        box[boxIndex] = temp_box;
        temp_parent.child[temp_box.child_index] = 0;
        box[temp_parent.boxid - 1] = temp_parent;
        pthread_mutex_unlock(&lock_box);
    }
    else if (temp_box.n <= S)
    {         
        // cube is a leaf
        pthread_mutex_lock(&lock_leaf);
        leaf_counter++;
        leaf = (Box*)realloc(leaf, leaf_counter * sizeof(Box));
        leaf[leaf_counter - 1] = temp_box;
        // update total number of points (just to check)
        num_points += temp_box.n;
        pthread_mutex_unlock(&lock_leaf);
   
        pthread_mutex_lock(&lock_box);
        box[boxIndex] = temp_box;
        pthread_mutex_unlock(&lock_box);
    } 
    else
    { 
        pthread_mutex_lock(&lock_box);
        box[boxIndex] = temp_box;
        pthread_mutex_unlock(&lock_box);
        // create 8 subcubes
        cubeDivision(&temp_box);
    }
    return NULL;
}


void* cubeDivision(void *arg)
{
    Box cube = *(Box *)arg;
    int i, j, pos[8];
    
    pthread_mutex_lock(&lock_box);
    // allocate memory for 8 more (sub)cubes
    box = (Box*)realloc(box,(8 + box_counter) * sizeof(Box));

    // initialize subcubes (children)
    for(i = 0; i < 8; i++)
    {
        box_counter++;
        
        box[box_counter - 1].level = cube.level + 1;
        box[box_counter - 1].boxid = box_counter;
        box[box_counter - 1].parent = cube.boxid;
        box[box_counter - 1].length = cube.length / 2;
        box[box_counter - 1].n = 0;
        box[box_counter - 1].child_index = i;
        
        // update parent with his child
        box[cube.boxid - 1].child[i] = box_counter;
        
        // initialize colleagues
        box[box_counter - 1].colleague_counter = 0;
        for (j = 0; j < 26; j++)
        {
            box[box_counter - 1].colleague[j] = 0;
        }
    }   
    
    cube.temp_counter = box_counter; 
    
    /* Set subcubes centers*/
    // Left - Front - Down
    box[box_counter - 8].center[0] = cube.center[0] - cube.length / 4;
	box[box_counter - 8].center[1] = cube.center[1] - cube.length / 4;
	box[box_counter - 8].center[2] = cube.center[2] - cube.length / 4;

	// Left - Front - Up
	box[box_counter - 7].center[0] = cube.center[0] - cube.length / 4;
	box[box_counter - 7].center[1] = cube.center[1] - cube.length / 4;
	box[box_counter - 7].center[2] = cube.center[2] + cube.length / 4;

	// Left - Back - Down
	box[box_counter - 6].center[0] = cube.center[0] - cube.length / 4;
	box[box_counter - 6].center[1] = cube.center[1] + cube.length / 4;
	box[box_counter - 6].center[2] = cube.center[2] - cube.length / 4;

	// Left - Back - Up
	box[box_counter - 5].center[0] = cube.center[0] - cube.length / 4;
	box[box_counter - 5].center[1] = cube.center[1] + cube.length / 4;
	box[box_counter - 5].center[2] = cube.center[2] + cube.length / 4;

	// Right - Front - Down
	box[box_counter - 4].center[0] = cube.center[0] + cube.length / 4;
	box[box_counter - 4].center[1] = cube.center[1] - cube.length / 4;
	box[box_counter - 4].center[2] = cube.center[2] - cube.length / 4;

	// Right - Front - Up
	box[box_counter - 3].center[0] = cube.center[0] + cube.length / 4;
	box[box_counter - 3].center[1] = cube.center[1] - cube.length / 4;
	box[box_counter - 3].center[2] = cube.center[2] + cube.length / 4;

	// Right - Back - Down
	box[box_counter - 2].center[0] = cube.center[0] + cube.length / 4;
	box[box_counter - 2].center[1] = cube.center[1] + cube.length / 4;
	box[box_counter - 2].center[2] = cube.center[2] - cube.length / 4;

	// Right - Back - Up
	box[box_counter - 1].center[0] = cube.center[0] + cube.length / 4;
	box[box_counter - 1].center[1] = cube.center[1] + cube.length / 4;
	box[box_counter - 1].center[2] = cube.center[2] + cube.length / 4;
    
    // check if we have new max level
    if (cube.level + 1 > num_levels)
    {
        num_levels = cube.level + 1;
    } 
    pthread_mutex_unlock(&lock_box);

    for (i = 0; i < 8; i++)
    {
		pos[i] = cube.temp_counter - i - 1;
	}
    if (running_threads < num_threads)
    { 
        // create 8 new threads
        pthread_mutex_lock(&lock_thread);
        running_threads += 8;
        pthread_mutex_unlock(&lock_thread);
        
        pthread_t t1, t2, t3, t4, t5, t6, t7, t8;

        pthread_create(&t1, NULL, &cubeCheck, (void *)&pos[7]);
        pthread_create(&t2, NULL, &cubeCheck, (void *)&pos[6]);
        pthread_create(&t3, NULL, &cubeCheck, (void *)&pos[5]);
        pthread_create(&t4, NULL, &cubeCheck, (void *)&pos[4]);
        pthread_create(&t5, NULL, &cubeCheck, (void *)&pos[3]);
        pthread_create(&t6, NULL, &cubeCheck, (void *)&pos[2]);
        pthread_create(&t7, NULL, &cubeCheck, (void *)&pos[1]);
        pthread_create(&t8, NULL, &cubeCheck, (void *)&pos[0]);

        // block the calling thread until the specified threads terminate
        pthread_join(t1, NULL);
        pthread_join(t2, NULL);
        pthread_join(t3, NULL);
        pthread_join(t4, NULL);
        pthread_join(t5, NULL);
        pthread_join(t6, NULL);
        pthread_join(t7, NULL);
        pthread_join(t8, NULL);

         pthread_mutex_lock(&lock_thread);
         // reduce number of running threads when these 8 terminate
         running_threads -= 8;
         pthread_mutex_unlock(&lock_thread);

         free(cube.points);
    }
    else 
    {
        // we can't create more threads
        for (i = 7; i >= 0; i--)
        {
            cubeCheck(&pos[i]);
        }
    }
    return NULL;
}


void searchNeighbours()
{
    int level, i, j,m,parent_id,child_id,colleague_id,colleague_index;
    double dist0,dist1,dist2;

    /* find colleagues searching level by level */
    for (level = 0; level < num_levels + 1; level++)
    {
        // search in all boxes
        for (i = 1; i < box_counter; i++)
        {
            if (box[i].level == level)
            {
                parent_id = box[i].parent;
                if (parent_id != 0)
                {
                    for (j = 0; j < 8; j++)
                    {
                        child_id = box[parent_id - 1].child[j];
                        if (child_id != 0)
                        {
                            if (box[i].boxid != box[child_id - 1].boxid)
                            {
                                // all "brothers" are colleagues (we can ignore the distance)
                                // we found a colleague!
                                box[i].colleague[box[i].colleague_counter++] = box[child_id - 1].boxid;
                            }
                        }
                    }
            
                    for (j = 0; j < 26; j++)
                    {
                        colleague_id = box[parent_id - 1].colleague[j]; 
                        // check if parent's colleague has children (if it's not empty or a leaf one) 
                        if (colleague_id != 0)
                        {
                            if (box[colleague_id - 1].n > S)
                            {    
                                for (m = 0; m < 8; m++)
                                {
                                    child_id = box[colleague_id - 1].child[m];
                                    if (child_id != 0)
                                    {
                                        if (box[i].boxid != box[child_id - 1].boxid)
                                        {
                                            // calculate distances
                                            dist0 = box[child_id - 1].center[0] - box[i].center[0];
                                            dist1 = box[child_id - 1].center[1] - box[i].center[1];
                                            dist2 = box[child_id - 1].center[2] - box[i].center[2]; 
                                        
                                            //check if distance is <=root(3)*length (= for the case when we have one common point only)
                                            if(sqrt(dist0 * dist0 + dist1 * dist1 + dist2 * dist2) <= sqrt(3) * box[i].length)
                                            {    
                                                colleague_index = box[i].colleague_counter;
                                                box[i].colleague[colleague_index] = box[child_id - 1].boxid;
                                                box[i].colleague_counter++;
                                            }
                                        }
                                    }
                    
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}  


void checkB()
{
    int i, j, same_counter = 0;
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++)
        {
            if ((B[i][0] == A[j][0]) && (B[i][1] == A[j][1]) && (B[i][2] == A[j][2]))
            {
                same_counter++;
            }
        }
    }
    if (same_counter == N)
    {
        printf("All points of B are also points of A\n");
    }
    else
    {
        printf("Error with points of B\n");
    }
}


void checkBoundaries()
{
    int i, j, points_counter = 0;
    double x, y, z;
	
    for (i = 0; i < leaf_counter; i++)
    {
        for (j = 0; j < leaf[i].n; j++)
        {
            x = fabs(leaf[i].center[0] - A[leaf[i].points[j]][0]);
            y = fabs(leaf[i].center[1] - A[leaf[i].points[j]][1]);
            z = fabs(leaf[i].center[2] - A[leaf[i].points[j]][2]);
            if (x < leaf[i].length / 2 && y < leaf[i].length / 2 && z < leaf[i].length / 2)
            {
                points_counter++;
            }
        }
    }
    if (points_counter == N)
    {
        printf("All points of leafs meet boundaries of subcubes\n");
    }
    else
    {
        printf("Error with points of leafs\n");
    }
}


//write array A to file
void AFile()
{
    remove("A.txt");
    
    FILE *A_file;
    int i;
    
    A_file = fopen("A.txt","wt");
    
    for (i = 0; i < N; i++)
    {
        fprintf(A_file,"%f,%f,%f\n",A[i][0],A[i][1],A[i][2]); fflush(A_file);
    }
    
    fclose(A_file);
}


//write array B to file
void BFile()
{
    remove("B.txt");
    
    FILE *B_file;
    int i;
    
    B_file = fopen("B.txt","wt");
    
    for (i = 0; i < N; i++)
    {
        fprintf(B_file,"%f,%f,%f\n",B[i][0],B[i][1],B[i][2]); fflush(B_file);  
    }
    
    fclose(B_file);
}


//write neighbours to file
void neighboursFile()
{   
    remove("neighbours.txt");
    
    FILE *neighbours_file;
    int i,j;
    
    neighbours_file = fopen("neighbours.txt","wt");
    for (i = 0; i < box_counter; i++)
    {               
        if (box[i].boxid != 0)
        {
            fprintf(neighbours_file,"id: %8d   neighbours:",box[i].boxid); fflush(neighbours_file);
            for (j = 0; j < 26; j++)
            {
                if (box[i].colleague[j] != 0)
                {
                    fprintf(neighbours_file,"%8d",box[i].colleague[j]); fflush(neighbours_file);
                }
            }
            fprintf(neighbours_file,"\n"); fflush(neighbours_file);
        }   
    }
    
    fclose(neighbours_file);   
}