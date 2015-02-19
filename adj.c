/*
 * David Teitelbaum ENEE446 assignment 1 problem 2 submission.
 * 02/19/2015
 *
 * My algorithm using the Simulated Annealing method to search
 * for the most unique pairs. With this method, I can consistently
 * obtain 197 unique pairs @ 1M samples (or about 3 minutes).
 *
 * Run 'Make' from any Mac or Linux terminal to build.
 *
 * 197 pair sequence
 *
 *
 * int seq_mat[20][5] = {{12, 5, 13, 15, 18 },                                                                                                                                                                                                     
 *                    {6, 7, 17, 3, 18 },                                                                                                                                                                                                      
 *                    {7, 14, 4, 13, 1 },                                                                                                                                                                                                      
 *                    {1, 2, 15, 17, 20 },                                                                                                                                                                                                       
 *                    {20, 14, 0, 12, 3 },                                                                                                                                                                                                    
 *                    {19, 8, 3, 5, 1 },                                                                                                                                                                                                      
 *                     {6, 11, 16, 12, 1 },                                                                                                                                                                                                      
 *                    {16, 10, 3, 15, 4 },                                                                                                                                                                                                      
 *                    {15, 0, 19, 11, 7 },                                                                                                                                                                                                     
 *                    {10, 11, 5, 17, 14 },                                                                                                                                                                                                     
 *                    {9, 19, 17, 12, 4},                                                                                                                                                                                                      
 *                    {8, 0, 13, 17, 16},                                                                                                                                                                                                      
 *                    {10, 2, 12, 7, 8   },                                                                                                                                                                                                     
 *                    {7, 20, 9, 5, 16},                                                                                                                                                                                                    
 *                    {10, 6, 20, 19, 13 },                                                                                                                                                                                                     
 *                    {18, 20, 11, 4, 8},                                                                                                                                                                                                   
 *                    {0, 18, 9, 1, 10},                                                                                                                                                                                                      
 *                    {16, 2, 14, 6, 19 },                                                                                                                                                                                                      
 *                    {13, 3, 2, 9, 11},
 *                    {5, 2, 4, 6, 0 }};
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "glob.h"
#include "adj.h"

/*
 * Tunable params
 */
#define NUM_ITR 1000000.
#define FINAL_TEMP .000001
#define INIT_TEMP 2.75

#define WITH_CLAMP_OPT
#define CLAMP .0007


/*
 * Returns an acceptance probability
 * given a current e, neighboring e,
 * and temperature
 */
float
get_prob(int e, int e_p, float temp)
{
    float prob = 0;

    // invert energies
    float e_p_f = (210 - (e_p));
    float e_f = (210 - (e));     

    if (e_p_f < e_f) {
        prob = 1;
    } else {
        prob = exp(0 - ((e_p_f - e_f)/temp));
#ifdef WITH_CLAMP_OPT
        // slight opt. I'm clamping acceptance prob
        if (prob < CLAMP) { 
            prob = CLAMP;
        }
#endif
    }
    return prob;
}

/*
 * Returns a temp given an index in the main 
 * SA loop.
 *
 */
float
get_temp(float t0, float t_end, float i, int N)
{
    // This is tunable too. I'm picking the best
    // temp decay function I had...
    // 194 unique pairs @ 500k samples 
    float a = log(t0 - t_end) / (float) log(N);
    return  t0 - pow(i, a);

}

int 
contains(int *seq, int ind)
{
    int i;
    for (i = 0; i < NUM_K; i++) {
        if ( seq[i] == ind ) {
            return 1;
        }
    }
    return 0;
}

/*
 * Given a current state, get a random neighbor,
 * where a neighbor can be a new pair in an any
 * of the patients.
 *
 * We only move one step away from the current
 * state (e.g. one pair change)
 */
int **
get_neighbor(int **seq, int **adj)
{
    int **n = (int **)malloc(sizeof(int *) * NUM_PATIENTS);
    int i;

    for (i = 0; i < NUM_PATIENTS; i++) {
        n[i] = (int *)malloc(sizeof(int) * NUM_K);
        // copy errthing ova
        memcpy(n[i], seq[i], sizeof(int) *  NUM_K); 
    }    
    // get random patient
    int patient_ind = rand_limit(NUM_PATIENTS);
    
    // get random bit to flip in current sequence
    int rand_ind = rand_limit(NUM_K);
    int drug_ind = rand_ind;
    
    // get another random bit to flip for the new seq
    rand_ind = rand_limit(NUM_N - NUM_K); // 0 -> 6
    int rand_counter = 0;

    for (i = 0; i < NUM_DRUGS; i++) {
        if (adj[patient_ind][i] == 1) { // we have an entry in seq
            if (!contains(seq[patient_ind], i)) {
                if (rand_counter == rand_ind) {
                    n[patient_ind][drug_ind] = i;
                    return n;
                } else {
                    rand_counter++;
                }
            }
        }
    }
    return NULL;    
}

int **
get_adj_matrix()
{
    int i, j;
    int **adj = (int **) calloc(sizeof(int *) * NUM_PATIENTS, 1);
    
    for (i = 0; i < NUM_PATIENTS; i++) {
        adj[i] = (int *) calloc(sizeof(int) * NUM_DRUGS, 1);
    }       
    for (i = 0; i < NUM_PATIENTS; i++) {
        for (j = 0; j < NUM_DRUGS; j++) {
            adj[i][j] = _adj[i][j];
        }
    }
    return adj; 
}

/*
 * Prints adjacency matrix
 */
void
print_adj(int **adj)
{
    int i, j;
    for (i = 0; i < NUM_PATIENTS; i++) {
        printf("%d: ", i);
        for (j = 0; j < NUM_DRUGS; j++) {
            printf("%d ", adj[i][j]);            
        }
        printf("\n");
    }
}

void
print_seq(int **seq)
{
   int i, j;
   for (i = 0; i < NUM_PATIENTS; i++) {
       for (j = 0; j < NUM_K; j++) {
           printf("%d ", seq[i][j]);            
       }
       printf("\n");
   }
   printf("==========\n");
}

/*
 * Generic 'get N choose K set' func
 */
int itr;
void subset(int arr[], int data[], int start, int end, 
            int index, int r, int **r_seqs)
{
    int i, j;
    if (index == r) {
        for (j = 0; j < r; j++) {
            r_seqs[itr][j] = data[j];
        } itr++;
        return;
    }
    for (i = start; i <= end && end - i + 1 >= r - index; i++) {
        data[index] = arr[i];
        subset(arr, data, i+1, end, index+1, r, r_seqs);
    }
}

int **
get_row_seqs(int *row)
{
    int i, j=0;    
    int row_indices[NUM_N] = {0};
    int **r_seqs = (int **)malloc(sizeof(int *) * N_C_K);    
    int buf[NUM_N];

    for (i = 0; i < N_C_K; i++) {
        r_seqs[i] = (int *)malloc(sizeof(int *) * NUM_K);       
    }
    for (i = 0; i < NUM_DRUGS; i++) {
        if (row[i] == 1) {
            row_indices[j++] = i;
        }
    }
    /*
     * Recursive func that generates 11 choose 5 combinations 
     * into r_seqs
     */
    itr = 0;
    subset(row_indices, buf, 0, NUM_N - 1, 0, NUM_K, r_seqs);
    return r_seqs;
}

/*
 * Prepopulate all sequences for each patient
 */
int ***
prepop_seqs( int **adj )
{
    
    int ***seqs = (int ***)malloc(sizeof(int **) * NUM_PATIENTS);
    int i;

    for (i = 0; i < NUM_PATIENTS; i++) {
        seqs[i] = get_row_seqs(adj[i]);
    }        
    return seqs;
}

/*
 * Given a sequence, count the number of pairs
 */
int 
calc_pairs( int **seq )
{
    char *big_tbl = (char *)calloc(2097152, 1); // faster than hashing!
    int i, j, k;
    int pair_cnt = 0;
    int pair_first, pair_sec;    

    for (i = 0; i < NUM_PATIENTS; i++) { 
        for (j = 0; j < NUM_K - 1; j++) {            
            for (k = j + 1; k < NUM_K; k++) {
                unsigned int tbl_index_1 = 
                    ((1 << seq[i][j]) | (1 << seq[i][k])); 

                unsigned int tbl_index_2 = 
                    ((1 << seq[i][k]) | (1 << seq[i][j])); 

                if (big_tbl[tbl_index_1] == 0 && 
                    big_tbl[tbl_index_2] == 0) {
                    big_tbl[tbl_index_1] = 1;
                    big_tbl[tbl_index_2] = 1;
                    pair_cnt++;
                }
            }
        }
    }
    free(big_tbl);
    return pair_cnt;
}

/*
 * Run the simulated annealing algorithm. Given some random
 * starting point, pick a random neighbor and start moving
 * around the system. As the system cools (temp starts high
 * and gets lower), acceptance probability will go down as
 * well, meaning the SA algorithm will converge to a simple
 * greedy hill climber (always move to the next best spot).
 *
 * This algorithm is great because you get the benefits of
 * the greedy alg without getting stuck in local maxima.
 */
int **
simulate_anneal( int **adj )
{
    int i,j;

    long int last = 0;
    long int now = 0;
    long int num_samples = 0;
    int high_num_pairs = 0;
    int *seq[NUM_PATIENTS];

    /*
     * Pre-populate all of the permutations of each row
     * in a 3D array of [20][126][5]
     */
    int ***all_seqs = prepop_seqs(adj);
    
    int max_score = 0;
    
    for (i = 0; i < NUM_PATIENTS; i++) {
        int rand = rand_limit(N_C_K);
        seq[i] = all_seqs[i][rand];
    }    
    
    max_score = calc_pairs(seq);
    printf("Max score is %d\n", max_score);
    for(i = 0; i < NUM_ITR; i++) {
        float temp = get_temp(INIT_TEMP, FINAL_TEMP, i, NUM_ITR);
        float acceptance_prob = 0;
        int **neighbor = get_neighbor(seq, adj);
        int score = calc_pairs(neighbor);
        
        if (score > max_score) {
            print_seq(neighbor);
            max_score = score;
            printf("Max score is %d\n", max_score);
        }
        acceptance_prob = get_prob(calc_pairs(seq), score, temp);
        
        if (acceptance_prob > ((float)rand() / (float)RAND_MAX)) {                 
            for (j = 0; j < NUM_PATIENTS; j++)
                seq[j] = neighbor[j];
        } else {
            for(j = 0; j < NUM_PATIENTS; j++)
                free(neighbor[j]);
                free(neighbor);
        }
        printf("Acceptance prob: %f temp: %f\r", 
               acceptance_prob, temp);
    }
    
    return adj;
}
