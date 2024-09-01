//
//Here, we define a class Solver_voles
//

#ifndef _solver_h
#define _solver_h

#include "computation_data.h"

enum Interaction{Repulsion, Attraction, Alignment};

class Solver: public Computation_data
{
private:
    VECTOR U_plus_current; //The density of right-moving with size (Nx) at current step
    VECTOR U_plus_next; //The density of right-moving with size (Nx) next step

    VECTOR U_minus_current; //The density of left-moving with size (Nx) at current step
    VECTOR U_minus_next; //The density of left-moving with size (Nx) next step

    //VECTOR Total_density_current; //The total dendity as u_plus + u_minus with size (Nx) at current step

    VECTOR Lambda_plus; //Turning rate: lambda^+

    VECTOR Lambda_minus; //Turning rate: lambda^-

    VECTOR Source_term_plus; // The source term for u_plus

    VECTOR Source_term_minus; // The source term for u_minus

    //store the result in files
    //int time_scale = 128;
    //int space_scale = 8;
    
    int number_of_final_step = 220;
    MATRIX U_plus_stable; //Save density of right-moving some final steps with size (number_of_final_step, Nx)
    MATRIX U_minus_stable; //Save density of left-moving some final steps with size (number_of_final_step, Nx)

    // int three_quarter;
    // MATRIX U_plus_three_quarter;
    // MATRIX U_minus_three_quarter;

    bool parallel;
    int number_of_threads; //the number of CPU is used in parallel computing. if = 0 then using default 4 CPUs

    VECTOR U_plus_temp, U_minus_temp; 
    VECTOR U_plus_current_save, U_minus_current_save; 


    VECTOR Runge_Kutta_1p, Runge_Kutta_2p, Runge_Kutta_3p, Runge_Kutta_4p;
    VECTOR Runge_Kutta_1m, Runge_Kutta_2m, Runge_Kutta_3m, Runge_Kutta_4m;

    //const double Tolerance = pow(10, -14);
    bool satisfied_tolerance = false;

    bool non_convergence = true;
    double non_convergence_error_min = pow(10, -8);
    double non_convergence_stop_time = 33000;

    double cpu_time_hour = 0;

    double CFL = 0;

    VECTOR L1_error;
    MATRIX Total_desity_100;

    int final_time_input;
    int Nt_input;

public:
    /**
     * @brief Construct a new Solver object with input parameters file
     * @param filename The name of the input parameters file
     */
    Solver(string filename);
    ~Solver();

    /**
     * @brief this function solve our problem
     */
    void Solve();

    /**
     * @brief solve the problem using parallel computing.
     */
    void Solve_parallel();

    /**
     * @brief solve the problem using parallel computing.
     * @param number_threads the number of CPU is used
     */
    void Solve_parallel(int number_threads);

    
    //VECTOR get_U_plus();
    //VECTOR get_U_minus();
    VECTOR get_Total_density_final_time();

    int get_Nt();
    double get_dt();

    int get_Nx();
    double get_dx();

    int get_Nr();
    int get_Na();
    int get_Nal();

    double get_CPU_time();
    double get_initial_amplitude();

    void set_initial_amplitude(double new_initial_amplitude);
    void set_scheme(Scheme new_scheme);
    void set_slope_limiter(SlopeLimiter new_slope);
    void set_q_a(double new_q_a);
    void set_output_directory(string new_output_directory);
    void set_final_time(double new_final_time);
    void set_dx(double new_dx);
    void set_dt(double new_dt);
    void set_tolerance(double new_tolerance);

private:
    /**
     * @brief check the correctness of input data
     */
    void check_all_data();

    /**
     * @brief write our solution to files
     */
    void write_solution();

    /**
     * @brief solve the problem 
     */
    void Solve_detail();

    /**
     * @brief compute our scheme to find a solution which is called in Solve()
     */
    void compute_solution();


    /**
     * @brief compute the integral in the signals perceived from neighbors: y^{+-} at current step
     * @param interaction Repulsion or Attraction or Alignment
     * @param i the position of U
     * @return double the value of the integral
     */
    double compute_integral(Interaction interaction, int i);

    /**
     * @brief compute the turning rates: Lambda_plus and Lambda_minus at current step
     */
    void compute_turning_rate();

    /**
     * @brief compute thr turning rates: Lambda_plus and Lambda_minus at current step using parallel computing 
     */
    void compute_turning_rate_parallel();

    /**
     * @brief compute K_j(s)(.....) depend on the model at current step
     * @param interaction Repulsion or Attraction or Alignment
     * @param i the position of U
     * @param s variable of the integral
     * @return double the value of K_j(s)(.....)
     */
    double compute_interaction_kernel_extension(Interaction interaction, int i, int s);

    /**
     * @brief Get the value of U_plus_current at the position depend on the boundary condition
     * @param position The position from -N_x to 2N_x-1
     * @return double the value of U_plus_current at the position
     */
    double get_U_plus_current(int position);

    /**
     * @brief Get the value of U_minus_current at the position depend on the boundary condition
     * @param position The position from -N_x to 2N_x-1
     * @return double the value of U_minus_current at the position
     */
    double get_U_minus_current(int position);

    /**
     * @brief Get the value of Total_density_current at the position depend on the boundary condition
     * @param position The position from -N_x to 2N_x-1
     * @return double the value of Total_density_current at the position
     */
    //double get_total_density_current(int position);

    /**
     * @brief Get the value of turning rate Lambda_plus at the position depend on the boundary condition
     * @param position The position from -N_x to 2N_x-1
     * @return double the value of turning rate Lambda_plus
     */
    double get_Lambda_plus(int position);

    /**
     * @brief Get the value of turning rate Lambda_minus at the position depend on the boundary condition
     * @param position The position from -N_x to 2N_x-1
     * @return double the value of turning rate Lambda_minus
     */
    double get_Lambda_minus(int position);

    double get_Source_plus(int position);

    double get_Source_minus(int position);

    /**
     * @brief Create objects for solver before computations to ensure that the memory of computer is enough
     */
    void create_objects();

    /**
     * @brief delete all old output files in output folder
     */
    void delete_old_output();

    /**
     * @brief compute the discrete error of total density between the current step and the next step
     * @return double the discrete error of total density between the current step and the next step
     */
    double compute_discrete_error();

    /**
     * @brief use the Upwind scheme to solve our problem
     */
    void Upwind_scheme();


    /**
     * @brief use the Downwind scheme to solve our problem
     */
    void Downwind_scheme();

    /**
     * @brief the high resoluion method using MC limiter with non-source term
     */
    void High_resoluion_method_non_source_term();

    /**
     * @brief use the Fractional step method to solve our problem with high resoluion method using MC limiter and Runge Kutta 2 ODE-solver
     */
    void Fractional_step_method_limiter_Runge_Kutta_scheme();

    /**
     * @brief use the MacCormack scheme to solve our problem
     */
    void MacCormack_scheme();

    /**
     * @brief use the Quasi-Steady Wave-Propagation Algorithm with slope limiter to solve our problem
     */
    void Quasi_Steady_Algorithm_limiter_scheme();

    /**
     * @brief use the Quasi-Steady Wave-Propagation Algorithm (non limiter) to solve our problem
     */
    void Quasi_Steady_Algorithm_scheme();

    /**
     * @brief use the Quasi-Steady Wave-Propagation Algorithm with centered slope to solve our problem
     */
    void Quasi_Steady_Algorithm_Center_scheme();

    /**
     * @brief use the Quasi-Steady Wave-Propagation Algorithm with BW slope to solve our problem
     */
    void Quasi_Steady_Algorithm_Beam_Warming_scheme();

    /**
     * @brief use the Quasi-Steady Wave-Propagation Algorithm with LW slope to solve our problem
     */
    void Quasi_Steady_Algorithm_Lax_Wendroof_scheme();

    /**
     * @brief use the F-Wave Algorithm to solve our problem
     */
    void F_Wave_Algorithm_scheme();

    /**
     * @brief check the tolerance condition at the time step n
     * @param n the time step n
     */
    void check_tolerance_condition(int n);

    /**
     * @brief compute two source term of our problem
     */
    void compute_source_term();

    /**
     * @brief copy value from U_next to U_current for preparing the next step
     */
    void copy_value_next_to_current();

    /**
     * @brief copy value from U_save to U_current
     */
    void copy_value_save_to_current();

    /**
     * @brief copy value from U_current to U_save
     */
    void copy_value_current_to_save();


    /**
     * @brief Set the initial data
     */
    void set_initial_data();

    /**
     * @brief Runge Kutta method to solve ODE with source term = scale_source*source_term
     * @param order The order: only 2 or 4
     */
    void Runge_Kutta_method(int order);

    /**
     * @brief compute the MC limiter at position i
     * @param i The position
     * @return double The value of MC limiter at position i
     */
    double limiter_plus(int i);

    double limiter_minus(int i);

    double limiter_plus_quasi_steady(int i);

    double limiter_minus_quasi_steady(int i);

    double compute_slope(double slope_center, double slope_up, double slope_down, string plus_or_minus);

    void write_input_computation_data();

    void after_propagation(int n);
};



#endif