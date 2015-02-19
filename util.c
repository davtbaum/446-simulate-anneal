#include <stdio.h>
#include <stdlib.h>
#include <time.h>

/*
 * Generate a number from 0 to max, written as a replacement of
 * rand() % N, which only provides uniform distribution if N
 * is a power of 2. 
 *
 * Resource: http://stackoverflow.com/questions/822323
 */
int 
rand_limit(int n) 
{

    static int seeded = 0;
    if (!seeded) {
        srand(time(NULL));
        seeded = 1;
    }

    if ((n - 1) == RAND_MAX) {
            return rand();
    } else {    
        // Chop off all of the values that would cause skew...
        long end = RAND_MAX / n; // truncate skew
        end *= n;
    
        // ... and ignore results from rand() that fall above that 
        // limit. (Worst case the loop condition should succeed 50% 
        // of the time, so we can expect to bail out of this loop 
        // pretty quickly.)
        int r;
        while ((r = rand()) >= end);       
        return r % n;
    }
}
