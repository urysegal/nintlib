#include "nintlib.h"
#include <stdio.h>
#include <math.h>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/geometries.hpp>

#include "libslater/libslater.h"


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
calculate_one_value( int dim_num, double x[], double shift )
{   
    x[6] = x[6] - shift;
    
    Quantum_Numbers quantum_numbers = {1,0,0};

//    STO_Basis_Function_Info inf1( 0.252, quantum_numbers);
//    STO_Basis_Function_Info inf2( 0.952, quantum_numbers);

//    STO_Basis_Function f1(inf1, {0, 0, -0.14142136});
//    STO_Basis_Function f2(inf2, {0.70710678, 0, 0.56568542});
//    STO_Basis_Function f3(inf1, {0.50710678, 0, 0.86568542});
//    STO_Basis_Function f4(inf2, {0.50710678, 0, 0.86568542});

    STO_Basis_Function_Info inf1( 1.625, quantum_numbers);
    STO_Basis_Function_Info inf2( 1.625, quantum_numbers);

    STO_Basis_Function f1(inf1, {0, 0, 0});
    STO_Basis_Function f2(inf2, {0, 0, 0});
    STO_Basis_Function f3(inf1, {0, 0, 0});
    STO_Basis_Function f4(inf2, {0, 0, 0});

    center_t r1 = {x[0], x[1], x[2]};
    center_t r2 = {x[3], x[4], x[5]};

    return f1.evaluate_conjugate(r1) * f3.evaluate_conjugate(r2) *
           f2.evaluate(r1) * f4.evaluate(r2) *
            (1.0/sqrt(std::numbers::pi))*exp(-x[6]*x[6]*distance_squared(r1,r2));
}


const double range_max = 10;

int
main()
{
    int eval_num = 0;
    int dim_num = 7;
    
    double a[dim_num] = {0};
    double b[dim_num] = {0};

    for ( int i = 0 ; i < dim_num ; ++i ) {
        a[i] = -range_max;
        b[i] = range_max;
    }
	
    /* p5_nd */
    
    auto res = p5_nd(calculate_one_value, 7, a,b, &eval_num);

    printf("Result p5_nd: %15.15f+%15.15fi, %d points\n", res.real(), res.imag(), eval_num);
    
    /* Box_nd*/
# define ORDER 5
    int j;
    double wtab[ORDER] = {
    0.236926885056189087514264040720,
    0.478628670499366468041291514836,
    0.568888888888888888888888888889,
    0.478628670499366468041291514836,
    0.236926885056189087514264040720 };
  double wtab2[ORDER];
  double xtab[ORDER] = {
    -0.906179845938663992797626878299,
    -0.538469310105683091036314420700,
     0.0,
     0.538469310105683091036314420700,
     0.906179845938663992797626878299 };
  double xtab2[ORDER];
 
 /* adjust quadrature from [-1,1] to [a,b] */
 for ( j = 0; j < ORDER; j++ )
  {
    xtab2[j] = ( xtab[j]*a[0] );
  }
  for ( j = 0; j < ORDER; j++ )
  {
    wtab2[j] = a[0] * wtab[j];
  } 
   auto res2 = box_nd ( calculate_one_value, 7, a, ORDER, xtab2, wtab2, &eval_num );
   printf("Result box_nd: %15.15f+%15.15fi, %d points\n", res2.real(), res2.imag(), eval_num);
# undef ORDER
   
   /* Romberg */
   int ind;
   double tol = 0.001;
   int it_max = 3;
   int sub_num[dim_num] = {0};
   for ( int i = 0 ; i < dim_num ; ++i ) {
        sub_num[i] = 5;
    }

  auto res3 = romberg_nd ( calculate_one_value, a, b, dim_num, sub_num, it_max, tol,
    &ind, &eval_num );
printf("Result romberg_nd: %15.15f+%15.15fi, %d points\n", res3.real(), res3.imag(), eval_num);
    return 0;
}
