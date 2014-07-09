#include <pthread.h>

/******************** Points Struct ********************/

typedef struct
{
    int level, boxid, parent, child[8], n, start, colleague[26];
    int temp_counter, colleague_counter, child_index, *points;
    double center[3], length;
} Box;

// The struct timeval structure represents an elapsed time
struct timeval startwtime, endwtime;