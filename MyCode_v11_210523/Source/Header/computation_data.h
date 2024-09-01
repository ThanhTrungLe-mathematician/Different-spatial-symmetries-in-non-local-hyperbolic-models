//
//Here, we define a class Computation_Data which contains all input data for computation
//

#ifndef _computation_data_h
#define _computation_data_h

#include "utilities.h"
#include "input_parameter.h"

//enum Model{M1, M2, M3, M4, M5};
//enum Scheme{Upwind, MacCormack, Fractional_Step_Method, Quasi_Steady_Algorithm, F_Wave_Algorithm};
//enum SlopeLimiter{Non_limiter, Centered_slope, Upwind_slope, Downwind_slope, Beam_Warming, Lax_Wendroof, Minmod, Superbee, MC_limiter};
//enum BoundaryCondition{DBC, NBC, PBC};

class Computation_data
{
protected:
    string input_filename; //the name of input parameters file
    Input_parameter * Input; //contains all information from the input parameter file
protected:
    //space domain
    double left_boundary;     // the left boundary of the domain
    double right_boundary;    // the right boundary of the domain
    double dx;      //the step

    //time
    double final_time;  // the final time
    double dt;          // the time step

    double tolerance; //the tolerance of two adjacent steps

    double gamma; //the speed

    double lambda_1; //represent a base-line turning rate

    double lambda_2; //represent a bias turning rate

    double q_r; //Magnitude of repulsion

    double q_a; //Magnitude of attraction

    double q_al; //Magnitude of alignment

    //Fixed values
    const double y_0 = 2.0; //shift of the turning function

    const double s_a = 1.0; //the spatial region for attractive

    const double s_al = 0.5; //the spatial region for alignment

    const double s_r = 0.25; //the spatial region for repulsive

    const double m_a = s_a/8.0; //Width of attraction kernel

    const double m_al = s_al/8.0; //Width of alignment kernel

    const double m_r = s_r/8.0; //Width of repulsion kernel

    //computated values
    int Nt; //the number of time step
    int Nx; //the number of space step
    int Na; //the number of the step for 0..2s_a
    int Nr; //the number of the step for 0..2s_r
    int Nal; //the number of the step for 0..2s_al

    string output_directory; 

    double initial_amplitude; //the amplitude of the initial condition

    Model model; //the model is used such as M1, M2, M3, M4, M5.
    string Model_string[5] = {"M1", "M2", "M3", "M4", "M5"};

    Scheme scheme; //the scheme is used such as Upwind, MacCormack
    string Scheme_string[5] = {"Upwind", "MacCormack", "Fractional_Step_Method", "Quasi_Steady_Algorithm", "F_Wave_Algorithm"};

    SlopeLimiter slope_limiter; //the slope limiter is used
    string Slope_string[9] = {"Non_limiter", "Centered_slope", "Upwind_slope", "Downwind_slope", "Beam_Warming", "Lax_Wendroof", "Minmod", "Superbee", "MC_limiter"};

    BoundaryCondition boundary_condition; //Dirichlet (DBC), Neumann (NBC), Periodic(PBC)
    string Boundary_string[3] = {"Dirichlet", "Neumann", "Periodic"};


public:
    Computation_data(string filename);
    ~Computation_data();

    /**
     * @brief read input data from files
     * @param filename The name of input parameters file
     */
    void read_computation_data();
    
protected:
    /**
     * @brief write the input data before computation
     */
    void write_computation_data();

    /**
     * @brief check the correctness of computation data
     */
    void check_computation_data();
};

#endif
