#include "nintlib.h"
#include <stdio.h>

double
calculate_one_value( int dim_num, double x[] )
{
    double res;
    for ( int i = 0 ; i < dim_num ; ++i ) {
        res += x[i];
    }
    return res;
}

int
main()
{
    int eval_num = 0;
    double a[7] = {0};
    double b[7] = {1,1,1,1,1,1,1};

    double res = p5_nd(calculate_one_value, 7, a,b, &eval_num);

    printf("Result: %f, %d points\n", res, eval_num);
    return 0;
}
