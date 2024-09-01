//
//  utilities.h
//
//
//  This file contains all the includes , types, and all global constants and functions
//



#ifndef _utilities_h
#define _utilities_h

#define _USE_MATH_DEFINES

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <iomanip>
#include <cstdlib>
#include <sstream>
#include <limits>
#include <ctime> 
#include <algorithm>
#include <queue>
#include <exception>
#include <omp.h>
#include <stdio.h>

using namespace std;


#define EPSILON 1.e-14

#define PRECISION 15

#define INPUT_PATH "../Input/"

#define OUTPUT_PATH "../Output/"

typedef vector<double> VECTOR;
typedef vector<VECTOR> MATRIX;
typedef vector<MATRIX> TENSOR3D;
typedef vector<TENSOR3D> TENSOR4D;



enum Model{M1, M2, M3, M4, M5};
enum Scheme{Upwind, MacCormack, Fractional_Step_Method, Quasi_Steady_Algorithm, F_Wave_Algorithm};
enum SlopeLimiter{Non_limiter, Centered_slope, Upwind_slope, Downwind_slope, Beam_Warming, Lax_Wendroof, Minmod, Superbee, MC_limiter};
enum BoundaryCondition{DBC, NBC, PBC};


//global functions

void Randomize();

/**
 * @brief random number in [0,1)
 * @return double a random number in [0,1)
 */
double Rand();

/**
 * @brief the number of elements of a vector from start_point to end_point with the step
 * @param start_point The start point
 * @param end_point The end point
 * @param step The step
 * @return int the number of elements of the vector
 */
int Get_Length(double start_point, double end_point, double step);

/**
 * @brief create a vector from start point to end point with step
 * @param aVector The vector is created
 * @param start_point The start point
 * @param end_point The end point
 * @param step The step
 */
void Create_Vector(VECTOR &aVector, double start_point, double end_point, double step);

/**
 * @brief resize of a vector
 * @param aVector The vector is resized
 * @param Nx The new size
 */
void Resize(VECTOR &aVector, int Nx);

/**
 * @brief resize a matrix
 * @param aMatrix The matrix is resized
 * @param Nx The first parameter of new size
 * @param Ny The second parameter of new size
 */
void Resize(MATRIX &aMatrix, int Nx, int Ny);

/**
 * @brief resize a tensor3D
 * @param aTensor3D The tensor3D is resized
 * @param Nt The first parameter of new size
 * @param Nx The second parameter of new size
 * @param Ny The third parameter of new size
 */
void Resize(TENSOR3D &aTensor3D, int Nt, int Nx, int Ny);

/**
 * @brief resize a Tensor4D
 * @param aTensor4D The tensor4D is resized
 * @param Nt The first parameter of new size
 * @param Na The second parameter of new size
 * @param Nx The third parameter of new size
 * @param Ny The fourth parameter of new size
 */
void Resize(TENSOR4D &aTensor4D, int Nt, int Na, int Nx, int Ny); 

/**
 * @brief clean a vector
 * @param aVector The vector is cleaned
 */
void Destroy(VECTOR &aVector);

/**
 * @brief clean a matrix
 * @param aMatrix The matrix is cleaned
 */
void Destroy(MATRIX &aMatrix);

/**
 * @brief clean a tensor3D
 * @param aTensor3D The tensor3D is cleaned
 */
void Destroy(TENSOR3D &aTensor3D);

/**
 * @brief clean a tensor 4D
 * @param aTensor4D The tensor4D is cleaned
 */
void Destroy(TENSOR4D &aTensor4D);

/**
 * @brief print a vector
 * @param aVector The vector is printed
 */
void Print(VECTOR aVector);

/**
 * @brief print a matrix
 * @param aMatrix The matrix is printed
 */
void Print(MATRIX aMatrix);

/**
 * @brief print a tensor3D
 * @param aTensor3D The tensor3D is printed
 */
void Print(TENSOR3D aTensor3D);

/**
 * @brief write a vector to a file
 * @param aVector The vector is written
 * @param filename The name of the file
 */
void Write_Data(VECTOR &aVector, string filename);

/**
 * @brief write sum of two vectors to a file
 * @param aVector the first vector
 * @param bVector the second vector
 * @param filename the name of the file
 */
void Write_Data(VECTOR &aVector, VECTOR &bVector, string filename);

/**
 * @brief write a matrix to a file
 * @param aMatrix The matrix is written
 * @param filename The name of the file
 */
void Write_Data(MATRIX &aMatrix, string filename);

/**
 * @brief write sum of two matrixs to a file
 * @param aMatrix the first matrix
 * @param bMatrix the second matrix
 * @param filename the name of the file
 */
void Write_Data(MATRIX &aMatrix, MATRIX &bMatrix, string filename);

/**
 * @brief write a matrix to a file with a number of final step
 * @param aMatrix The matrix is written
 * @param filename The name of the file
 * @param number_of_final_step  The number of final step in x
 */
void Write_Data(MATRIX &aMatrix, string filename, int number_of_final_step);

/**
 * @brief write sum of two matrixs to a file with a number of final step in x
 * @param aMatrix the first matrix
 * @param bMatrix the second matrix
 * @param filename the name of the file
 * @param number_of_final_step the number of final step in x
 */
void Write_Data(MATRIX &aMatrix, MATRIX &bMatrix, string filename, int number_of_final_step);

/**
 * @brief write a matrix to a file with scales
 * @param aMatrix The matrix is written
 * @param filename The name of the file
 * @param scale_x The scale in x
 * @param scale_y The scale in y
 */
void Write_Data(MATRIX &aMatrix, string filename, int scale_x, int scale_y);

/**
 * @brief write sum of two matrixs to a file with scales
 * @param aMatrix the first matrix
 * @param bMatrix the second matrix
 * @param filename the name of the file
 * @param scale_x the scale in x
 * @param scale_y the scale in y
 */
void Write_Data(MATRIX &aMatrix, MATRIX &bMatrix, string filename, int scale_x, int scale_y);

/**
 * @brief write a tensor3D to a file
 * @param aTensor3D The tensor3D is written
 * @param filename The name of the file
 */
void Write_Data(TENSOR3D &aTensor3D, string filename);

/**
 * @brief read a vector from a file
 * @param aVector The vector is read
 * @param filename The name of the file
 */
void Read_Data(VECTOR &aVector, string filename);

/**
 * @brief read a matrix from a file
 * @param aMatrix The matrix is read
 * @param filename The name of the file
 */
void Read_Data(MATRIX &aMatrix, string filename);

/**
 * @brief read a tensor3D from a file
 * @param aTensor3D The tensor3D is read
 * @param filename The name of the file
 */
void Read_Data(TENSOR3D &aTensor3D, string filename);

/**
 * @brief find the maximum value of a vector
 * @param aVector The vector is use to find
 * @return double The maximum value of the vector
 */
double Get_Max_Value(VECTOR &aVector);

/**
 * @brief find the maximum value of a matrix
 * @param aMatrix The matrix is use to find
 * @return double The maximum value of the matrix
 */
double Get_Max_Value(MATRIX &aMatrix);

/**
 * @brief find the maximum value of a tensor3D
 * @param aTensor3D The tensor3D is use to find
 * @return double The maximum value of the tensor3D
 */
double Get_Max_Value(TENSOR3D &aTensor3D);

/**
 * @brief find the maximum value of a tensor4D
 * @param aTensor4D The tensor4D is use to find
 * @return double The maximum value of the tensor4D
 */
double Get_Max_Value(TENSOR4D &aTensor4D);


/**
 * @brief create a Vector depends on the function f(x)
 * @param aVector The vector is created
 * @param function The function which defines the value of f(x)
 * @param V_x The Vector constains all value of the variable x
 */
void Create_Vector(VECTOR &aVector, double (*function)(double), VECTOR V_x);


/**
 * @brief create a Matrix depends on the function f(x,y)
 * @param aMtrix The Matrix is created
 * @param function The function which defines the value of f(x,y)
 * @param V_x The Vector constains all value of the variable x
 * @param V_y The Vector constains all value of the variable y
 */
void Create_Matrix(MATRIX &aMtrix, double (*function)(double, double), VECTOR V_x, VECTOR V_y);


/**
 * @brief create a Tensor3D depends on the function f(t,x,y)
 * @param aTensor3D The Tensor3D is created
 * @param function The function which defines the value of f(t,x,y)
 * @param V_t The Vector constains all value of the variable time
 * @param V_x The Vector constains all value of the variable x
 * @param V_y The Vector constains all value of the variable y
 */
void Create_Tensor3D(TENSOR3D &aTensor3D, double (*function)(double, double, double), VECTOR V_t, VECTOR V_x, VECTOR V_y);


/**
 * @brief create a Tensor4D depends on the function f(t,a,x,y)
 * @param aTensor4D The Tensor4D is created
 * @param function The function which defines the value of f(t,a,x,y)
 * @param V_t The Vector constains all value of the variable time
 * @param V_a The Vector constains all value of the variable age
 * @param V_x The Vector constains all value of the variable x
 * @param V_y The Vector constains all value of the variable y 
 */
void Create_Tensor4D(TENSOR4D &aTensor4D, double (*function)(double, double, double, double), VECTOR V_t, VECTOR V_a, VECTOR V_x, VECTOR V_y);

/**
 * @brief compute the L1 discrete error of numerical solution compared with reference solution
 * @param Numerical_sol  The numerical solution
 * @param Reference_sol  The reference solution
 * @return double The L1 discrete error
 */
double compute_L1_discrete_error(TENSOR3D &Numerical_sol , TENSOR3D  &Reference_sol);

/**
 * @brief compute the L1 discrete error of numerical solution compared with reference solution
 * @param Numerical_sol The numerical solution
 * @param Reference_sol The reference solution
 * @return double The L1 discrete error
 */
double compute_L1_discrete_error(MATRIX &Numerical_sol, MATRIX &Reference_sol);

/**
 * @brief compute the L1 discrete error of two vectors
 * @param First_vector The first vector
 * @param Second_vector The second vector
 * @return double the L1 discrete error of two vectors
 */
double compute_L1_discrete_error(VECTOR &First_vector, VECTOR &Second_vector);


/**
 * @brief compute the L1 discrete error of (A+B)-(C+D)
 * @param A_vector 
 * @param B_vector 
 * @param C_vector 
 * @param D_vector 
 * @return double the L1 discrete error
 */
double compute_L1_discrete_error(VECTOR &A_vector, VECTOR &B_vector, VECTOR &C_vector, VECTOR &D_vector);

double min_of(double a, double b);
double min_of(double a, double b, double c);
double max_of(double a, double b);
double max_of(double a, double b, double c);
double minmod(double a, double b);
double minmod(double a, double b, double c);
double maxmod(double a, double b);
double maxmod(double a, double b, double c);

double get_cpu_time_sec();
double get_cpu_time_hour();
string get_local_time();

#endif
