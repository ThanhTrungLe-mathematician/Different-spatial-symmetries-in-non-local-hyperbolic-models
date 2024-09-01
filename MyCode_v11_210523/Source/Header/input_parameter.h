//
//Here, we define a class Input_parameter which include all data from input files
//

#ifndef _input_parameter_h
#define _input_parameter_h

#include "utilities.h" 

class Input_parameter
{
private:
    //space domain
    double left_boundary;     // the left boundary of the domain
    double right_boundary;    // the right boundary of the domain
    double dx;      //the step

    //time
    double final_time;  // the final time
    double dt;          // the time step
   
    string output_directory;    //The name of the output directory

    double gamma; //the speed

    double lambda_1; //represent a base-line turning rate

    double lambda_2; //represent a bias turning rate

    double q_r; //Magnitude of repulsion

    double q_a; //Magnitude ofattraction

    double q_al; //Magnitude of alignment

    double initial_amplitude; //the amplitude of the initial condition

    string model; //the model is used such as M1, M2, M3, M4, M5.

    string scheme; //the scheme is used such as Upwind, MacCormack

    string slope_limiter;

    string boundary_condition; //Boundary condition: Dirichlet (DBC), Neumann (NBC), Periodic(PBC)

    double tolerance; //Tolerance of two adjacent steps (10^)

private:
    //use for parser text file
    int number_of_lines;
    vector<string> data; // tabular that contains the strings "description <variable>:: value"

    /**
     * @brief parser text file with syntax: "description <variable>:: value"
     * @param filename the name of file to parser 
     */
    void Parser_textfile(string filename);

    /**
     * @brief return the value of the variable
     * @param name_of_variable "<name of the variable>"
     * @return string The value of the variable
     */
    string Get_Value_of_variable(string name_of_variable);
public:

    /**
     * @brief Construct a new Input_parameter object from input file
     * @param filename The name of file which contains input parameters
     */
    Input_parameter(string filename);
    
    Input_parameter();
    ~Input_parameter();

    /**
     * @brief Get the left boundary object
     * @return double left_boundary
     */
    double get_left_boundary();

    /**
     * @brief Get the right boundary object
     * @return double right_boundary
     */
    double get_right_boundary();

    /**
     * @brief Get the dx object
     * @return double dx
     */
    double get_dx();

    /**
     * @brief Get the final time object
     * @return double final_time
     */
    double get_final_time();

    /**
     * @brief Get the dt object
     * @return double dt
     */
    double get_dt();

    /**
     * @brief Get the output directory object
     * @return string output_directory
     */
    string get_output_directory();

    /**
     * @brief Get the gamma object
     * @return double gamma
     */
    double get_gamma();

    /**
     * @brief Get the lambda 1 object
     * @return double lambda_1
     */
    double get_lambda_1();

    /**
     * @brief Get the lambda 2 object
     * @return double lambda_2
     */
    double get_lambda_2();

    /**
     * @brief Get the q_r object
     * @return double q_r
     */
    double get_q_r();

    /**
     * @brief Get the q_a object
     * @return double q_a
     */
    double get_q_a();

    /**
     * @brief Get the q_al object
     * @return double q_al
     */
    double get_q_al();

    /**
     * @brief Get the initial amplitude object
     * @return double initial_amplitude
     */
    double get_initial_amplitude();

    /**
     * @brief Get the model object
     * @return string model
     */
    string get_model();

    /**
     * @brief Get the scheme object
     * @return string scheme
     */
    string get_scheme();

    /**
     * @brief Get the slope limiter object
     * @return string the slope limiter
     */
    string get_slope_limiter();

    /**
     * @brief Get the boundary condition object
     * @return string boundary condition
     */
    string get_boundary_condition();

    /**
     * @brief Get the tolerance object
     * @return int the tolerance
     */
    double get_tolerance();
};

#endif

