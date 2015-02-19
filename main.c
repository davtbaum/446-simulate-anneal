#include <stdio.h>
#include "glob.h"

int 
main(void)
{
    
    int **adj_matrix;

    adj_matrix = get_adj_matrix();
    print_adj(adj_matrix);
    
    simulate_anneal(adj_matrix);
}
