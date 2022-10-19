#include "nintlib.h"
#include <stdio.h>
#include <math.h>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/geometries.hpp>
#include <boost/math/special_functions.hpp>

#include "/home/urysegal/libslater/libslater.h"


using namespace slater;
namespace bg = boost::geometry;

double distance_squared(const center_t &A, const center_t &B)
{
    bg::model::point<double, 3, bg::cs::cartesian> A_point(A[0], A[1], A[2]);
    bg::model::point<double, 3, bg::cs::cartesian> B_point(B[0], B[1], B[2]);
    auto res = bg::distance(A_point,B_point);
    return res*res;
}


std::complex<double>
calculate_one_value( int dim_num, double x[] )
{
    Quantum_Numbers quantum_numbers = {2,1,0};

    STO_Basis_Function_Info inf1( 0.252, quantum_numbers);
    STO_Basis_Function_Info inf2( 0.952, quantum_numbers);

    STO_Basis_Function f1(inf1, {0, 0, -0.14142136});
    STO_Basis_Function f2(inf2, {0.70710678, 0, 0.56568542});
    STO_Basis_Function f3(inf1, {0.50710678, 0, 0.86568542});
    STO_Basis_Function f4(inf2, {0.50710678, 0, 0.86568542});

    center_t r1 = {x[0], x[1], x[2]};
    center_t r2 = {x[3], x[4], x[5]};

    return f1.evaluate_conjugate(r1) * f3.evaluate_conjugate(r2) *
           f2.evaluate(r1) * f4.evaluate(r2) *
            exp(-x[6]*distance_squared(r1,r2))*(1.0/sqrt(x[6]));
}

int
main()
{


    int eval_num = 0;
    double a[7] = {0};
    double b[7] = {1,1,1,1,1,1,1};

    auto res = p5_nd(calculate_one_value, 7, a,b, &eval_num);

    printf("Result: %f+%fi, %d points\n", res.real(), res.imag(), eval_num);
    return 0;
}
