#include <stdio.h>
#include <stdlib.h>

#include <time.h>
#include <sys/time.h>

#include <mach/clock.h>
#include <mach/mach.h>
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

long int
get_nano(void)
{
    struct timespec ts;

    clock_serv_t cclock;
    mach_timespec_t mts;
    host_get_clock_service(mach_host_self(), CALENDAR_CLOCK, &cclock);
    clock_get_time(cclock, &mts);
    mach_port_deallocate(mach_task_self(), cclock);
    ts.tv_sec = mts.tv_sec;
    ts.tv_nsec = mts.tv_nsec;

    return ts.tv_nsec;
}
