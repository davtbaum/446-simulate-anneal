#define NUM_PATIENTS  20
#define NUM_DRUGS     21
#define NUM_N         11
#define NUM_K         5
#define N_C_K         462

int **get_adj_matrix();
void print_adj(int **adj);
void print_seq(int **seq);
int **simulate_anneal(int **adj);
int rand_limit(int n);
long int get_nano(void);
