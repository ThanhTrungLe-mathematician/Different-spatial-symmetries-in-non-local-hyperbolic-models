//
// This file contains all methods of class Computation_data
//

#include "computation_data.h"

Computation_data::Computation_data(string filename)
{
    this->input_filename = filename;
}

void Computation_data::read_computation_data()
{
    // read input parameters via class Input_parameter
    Input = new Input_parameter(input_filename);

    this->left_boundary = Input->get_left_boundary();
    this->right_boundary = Input->get_right_boundary();
    this->dx = Input->get_dx();

    this->final_time = Input->get_final_time();
    this->dt = Input->get_dt();
    this->tolerance = pow(10, Input->get_tolerance());

    this->gamma = Input->get_gamma();
    this->lambda_1 = Input->get_lambda_1();
    this->lambda_2 = Input->get_lambda_2();
    this->q_a = Input->get_q_a();
    this->q_al = Input->get_q_al();
    this->q_r = Input->get_q_r();

    this->output_directory = Input->get_output_directory();

    this->initial_amplitude = Input->get_initial_amplitude();

    // compute Nt, Nx, Na, Nr, Nal
    this->Nt = Get_Length(0, this->final_time, this->dt);
    this->Nx = Get_Length(this->left_boundary, this->right_boundary, this->dx);
    this->Na = Get_Length(0, 2 * s_a, dx);
    this->Nr = Get_Length(0, 2 * s_r, dx);
    this->Nal = Get_Length(0, 2 * s_al, dx);
    if(this->Na % 2 == 0)
    {
        this->Na += 1;
    }
    if(this->Nr % 2 == 0)
    {
        this->Nr += 1;
    }
    if(this->Nal % 2 == 0)
    {
        this->Nal += 1;
    }

    //read the model
    bool get_success = false;
    for(int i =0; i < Model_string->size(); i++){
        if(Model_string[i] == Input->get_model())
        {
            this->model = static_cast<Model>(i);
            get_success = true;
        }
    }
    if(get_success == false)
    {
        cerr << "The model " << Input->get_model() << " does not exist. Please check again!" << endl;
        exit(EXIT_FAILURE);
    }

    //read the schem
    get_success = false;
    for(int i =0; i < Scheme_string->size(); i++){
        if(Scheme_string[i] == Input->get_scheme())
        {
            this->scheme = static_cast<Scheme>(i);
            get_success = true;
        }
    }
    if(get_success == false)
    {
        cerr << "The scheme " << Input->get_scheme() << " does not exist. Please check again!" << endl;
        exit(EXIT_FAILURE);
    }

    //read the slope limiter
    get_success = false;
    for(int i =0; i < Slope_string->size(); i++){
        if(Slope_string[i] == Input->get_slope_limiter())
        {
            this->slope_limiter = static_cast<SlopeLimiter>(i);
            get_success = true;
        }
    }
    if(get_success == false)
    {
        cerr << "The slope limiter " << Input->get_slope_limiter() << " does not exist. Please check again!" << endl;
        exit(EXIT_FAILURE);
    }

    //read the boundary condition
    get_success = false;
    for(int i =0; i < Boundary_string->size(); i++){
        if(Boundary_string[i] == Input->get_boundary_condition())
        {
            this->boundary_condition = static_cast<BoundaryCondition>(i);
            get_success = true;
        }
    }
    if(get_success == false)
    {
        cerr << "The boundary condition " << Input->get_boundary_condition() << " does not exist. Please check again!" << endl;
        exit(EXIT_FAILURE);
    }
}

Computation_data::~Computation_data()
{
    if (Input != NULL)
    {
        delete Input;
        Input = NULL;
    }
}

void Computation_data::write_computation_data()
{

    // create the name of file
    string filename =this->output_directory + "/input_computation_data.txt";

    ofstream writefile(filename, ios::out);

    // if the file can not be open then stop
    if (!writefile)
    {
        cerr << "Impossible to open the " << filename.c_str() << " file!" << endl;
        exit(EXIT_FAILURE);
    }

    writefile.precision(PRECISION);

    time_t now = time(0);
    string local_time = ctime(&now);

    writefile << "### Start at " + local_time << endl;
    writefile << endl;

    writefile << "#####################################################################" << endl;
    writefile << "# Input data for computation" << endl;
    writefile << "#####################################################################" << endl;
    writefile << endl;

    writefile << "The left boundary of the domain <left_boundary>::" << this->left_boundary << endl;
    writefile << endl;

    writefile << "The right boundary of the domain <right_boundary>::" << this->right_boundary << endl;
    writefile << endl;

    writefile << "The step <dx>::" << this->dx << endl;
    writefile << endl;

    writefile << "Final time <final_time>::" << this->final_time << endl;
    writefile << endl;

    writefile << "Time step <dt>::" << this->dt << endl;
    writefile << endl;

    writefile << "Tolerance of two adjacent steps <tolerance>::" << this->tolerance << endl;
    writefile << endl;

    writefile << "The speed <gamma>::" << this->gamma << endl;
    writefile << endl;

    writefile << "Base line turning rate <lambda_1>::" << this->lambda_1 << endl;
    writefile << endl;

    writefile << "Bias turning rate <lambda_2>::" << this->lambda_2 << endl;
    writefile << endl;

    writefile << "Magnitude of repulsion <q_r>::" << this->q_r << endl;
    writefile << endl;

    writefile << "Magnitude of attraction <q_a>::" << this->q_a << endl;
    writefile << endl;

    writefile << "Magnitude of alignment <q_al>::" << this->q_al << endl;
    writefile << endl;

    writefile << "Shift of the turning function <y_0>::" << this->y_0 << endl;
    writefile << endl;

    writefile << "The spatial region for attractive <s_a>::" << this->s_a << endl;
    writefile << endl;

    writefile << "The spatial region for alignment <s_al>::" << this->s_al << endl;
    writefile << endl;

    writefile << "The spatial region for repulsive <s_r>::" << this->s_r << endl;
    writefile << endl;

    writefile << "Width of attraction kernel <m_a>::" << this->m_a << endl;
    writefile << endl;

    writefile << "Width of alignment kernel <m_al>::" << this->m_al << endl;
    writefile << endl;

    writefile << "Width of repulsion kernel <m_r>::" << this->m_r << endl;
    writefile << endl;

    writefile << "The amplitude of the initial condition <initial_amplitude>::" << this->initial_amplitude << endl;
    writefile << endl;

    writefile << "The model is used: 0=M1, 1=M2, 2=M3, 3=M4, 4=M5 <model>::" << this->model << endl;
    writefile << endl;

    writefile << "The scheme is used: 0=Upwind, 1=MacCormack <scheme>::" << this->scheme << endl;
    writefile << endl;

    writefile << "Boundary condition: 0=Dirichlet, 1=Neumann, 2=Periodic <boundary_condition>::" << this->boundary_condition << endl;
    writefile << endl;

    writefile.close();
}

void Computation_data::check_computation_data()
{

    // check the boundary on the x-axis
    if (this->right_boundary <= this->left_boundary)
    {
        cerr << "The left boundary of the domain: " << this->left_boundary << endl;
        cerr << "The right boundary of the domain: " << this->right_boundary << endl;
        cerr << "The right boundary must be greater than the left boundary!" << endl;
        exit(EXIT_FAILURE);
    }

    // check the final time
    if (this->final_time <= 0)
    {
        cerr << "The final time: " << this->final_time << endl;
        cerr << "The final time must be positive!" << endl;
        exit(EXIT_FAILURE);
    }

    // check the x-axis step
    if (this->dx <= 0)
    {
        cerr << "The step: " << this->dx << endl;
        cerr << "The step must be positive!" << endl;
        exit(EXIT_FAILURE);
    }

    // check the time step
    if (this->dt <= 0)
    {
        cerr << "The time step: " << this->dt << endl;
        cerr << "The time step must be positive!" << endl;
        exit(EXIT_FAILURE);
    }

    // check the number of the step
    if (this->Na % 2 == 0)
    {
        cerr<< "The number of the step for 0..2s_a: " <<this->Na << endl;
        cerr<< "The number of the step for 0..2s_a must be odd!" <<endl;
        exit(EXIT_FAILURE);
    }

    if (this->Nr % 2 == 0)
    {
        cerr<< "The number of the step for 0..2s_r: " <<this->Nr << endl;
        cerr<< "The number of the step for 0..2s_r must be odd!" <<endl;
        exit(EXIT_FAILURE);
    }

    if (this->Nal % 2 == 0)
    {
        cerr<< "The number of the step for 0..2s_al: " <<this->Nal << endl;
        cerr<< "The number of the step for 0..2s_al must be odd!" <<endl;
        exit(EXIT_FAILURE);
    }
    
}