//
// This file contains all methods of class Solver_voles
//

#include "solver.h"

Solver::Solver(string filenmae) : Computation_data(filenmae) {}

Solver::~Solver()
{
    Destroy(this->U_plus_current);
    Destroy(this->U_minus_current);
    // Destroy(this->Total_density_current);
    Destroy(this->U_plus_next);
    Destroy(this->U_minus_next);
    Destroy(this->Lambda_plus);
    Destroy(this->Lambda_minus);
    Destroy(this->U_plus_temp);
    Destroy(this->U_minus_temp);
    Destroy(this->U_plus_stable);
    Destroy(this->U_minus_stable);
    //Destroy(this->U_plus_three_quarter);
    //Destroy(this->U_minus_three_quarter);
    Destroy(this->Source_term_plus);
    Destroy(this->Source_term_minus);
    Destroy(this->Runge_Kutta_1p);
    Destroy(this->Runge_Kutta_2p);
    Destroy(this->Runge_Kutta_3p);
    Destroy(this->Runge_Kutta_4p);
    Destroy(this->Runge_Kutta_1m);
    Destroy(this->Runge_Kutta_2m);
    Destroy(this->Runge_Kutta_3m);
    Destroy(this->Runge_Kutta_4m);
    Destroy(this->U_minus_current_save);
    Destroy(this->U_plus_current_save);
    Destroy(this->L1_error);
    Destroy(this->Total_desity_100);
}

void Solver::set_initial_data()
{
    // set the initial condition

    double L = this->right_boundary - this->left_boundary; // length of domain
    // double L = 10;

    for (int i = 0; i < this->Nx; i++)
    {
        //U_plus_current[i] = 1.0 + 0.5 * this->initial_amplitude * Rand();
        //U_minus_current[i] = 1.0 + 0.5 * this->initial_amplitude * Rand();

        // U_plus_current[i] = 1.0 + 0.5 * this->initial_amplitude * abs(sin(0.2 *M_PI * i * dx));
        // U_minus_current[i] = 1.0 + 0.5 * this->initial_amplitude * abs(sin(0.2 *M_PI * i * dx));

        U_plus_current[i] = 1.0 + 0.25 * this->initial_amplitude * (1 + sin(2.0 / L * M_PI * i * dx));
        U_minus_current[i] = 1.0 + 0.25 * this->initial_amplitude * (1 + sin(2.0 / L * M_PI * i * dx));

        //U_plus_current[i] = 1.0 + 0.25 * this->initial_amplitude * (1 + sin(4.0 / L * M_PI * i * dx));
        //U_minus_current[i] = 1.0 + 0.25 * this->initial_amplitude * (1 + sin(4.0 / L * M_PI * i * dx));

        // U_plus_current[i] = 1.0 + 0.25 * this->initial_amplitude * sin(0.2 *M_PI * i * dx);
        // U_minus_current[i] = 1.0 + 0.25 * this->initial_amplitude * sin(0.2 *M_PI * i * dx);

        // U_plus_current[i] = 1.0 + 0.25 * this->initial_amplitude * sin(0.6 * M_PI * M_PI * i * dx);
        // U_minus_current[i] = 1.0 + 0.25 * this->initial_amplitude * sin(0.6 * M_PI * M_PI * i * dx);

        // Total_density_current[i] = U_plus_current[i] + U_minus_current[i];
    }
}

void Solver::compute_solution()
{

    this->set_initial_data();

    string filename_initial;
    filename_initial = this->output_directory + "/density_right_moving_initial.txt";
    Write_Data(this->U_plus_current, filename_initial.c_str());
    filename_initial = this->output_directory + "/density_left_moving_initial.txt";
    Write_Data(this->U_minus_current, filename_initial.c_str());
    filename_initial = this->output_directory + "/total_density_initial.txt";
    Write_Data(this->U_plus_current, this->U_minus_current, filename_initial.c_str());

    // compute the density
    if (this->scheme == Upwind)
    {
        // cout << "Solver by using " << Scheme_string[this->scheme] << endl;
        for (int n = 0; n < this->Nt - 1; n++)
        {
            this->Upwind_scheme();
            this->after_propagation(n);
        }
    }
    else if (this->scheme == MacCormack)
    {
        // cout << "Solver by using " << Scheme_string[this->scheme] << endl;
        for (int n = 0; n < this->Nt - 1; n++)
        {
            this->MacCormack_scheme();
            this->after_propagation(n);
        }
    }
    else if (this->scheme == Fractional_Step_Method)
    {

        // cout << "Solver by using " << Scheme_string[this->scheme] << " and " << Slope_string[this->slope_limiter] << endl;
        for (int n = 0; n < this->Nt - 1; n++)
        {
            this->Fractional_step_method_limiter_Runge_Kutta_scheme();
            this->after_propagation(n);
        }
    }
    else if (this->scheme == Quasi_Steady_Algorithm)
    {
        if (this->slope_limiter == Non_limiter)
        {
            for (int n = 0; n < this->Nt - 1; n++)
            {
                this->Quasi_Steady_Algorithm_scheme();
                this->after_propagation(n);
            }
        }
        else if (this->slope_limiter == Centered_slope)
        {
            for (int n = 0; n < this->Nt - 1; n++)
            {
                this->Quasi_Steady_Algorithm_Center_scheme();
                this->after_propagation(n);
            }
        }
        else if (this->slope_limiter == Beam_Warming)
        {
            for (int n = 0; n < this->Nt - 1; n++)
            {
                this->Quasi_Steady_Algorithm_Beam_Warming_scheme();
                this->after_propagation(n);
            }
        }
        else if (this->slope_limiter == Lax_Wendroof)
        {
            for (int n = 0; n < this->Nt - 1; n++)
            {
                this->Quasi_Steady_Algorithm_Lax_Wendroof_scheme();
                this->after_propagation(n);
            }
        }
        else
        {
            for (int n = 0; n < this->Nt - 1; n++)
            {
                this->Quasi_Steady_Algorithm_limiter_scheme();
                this->after_propagation(n);
            }
        }
    }
    else if (this->scheme == F_Wave_Algorithm)
    {
        // cout << "Solver by using " << Scheme_string[this->scheme] << endl;
        for (int n = 0; n < this->Nt - 1; n++)
        {
            this->F_Wave_Algorithm_scheme();
            this->after_propagation(n);
        }
    }
    else
    {
        cerr << "Please code for another scheme!" << endl;
        exit(EXIT_FAILURE);
    }
}

void Solver::after_propagation(int n)
{
    this->check_tolerance_condition(n);

    for (int i = 0; i < 100; i++)
    {
        if (n == this->Nt_input * i / 100)
        {
            // cout<<i<<endl;
            for (int j = 0; j < this->Nx; j++)
            {
                this->Total_desity_100[i][j] = this->U_plus_current[j] + this->U_minus_current[j];
            }
            break;
        }
    }

    this->copy_value_next_to_current();

    if (n % int(50000/this->dt) == 0)
    {
        this->write_solution();
    }

    if (n > this->Nt - 2 - this->number_of_final_step)
    {
        int time = n - (Nt - 2 - this->number_of_final_step) - 1;
        for (int i = 0; i < this->Nx; i++)
        {
            this->U_plus_stable[time][i] = this->U_plus_current[i];
            this->U_minus_stable[time][i] = this->U_minus_current[i];
        }
    }

    // if (n > three_quarter - this->number_of_final_step && n <= three_quarter)
    // {
    //     int time = n - (three_quarter - this->number_of_final_step) - 1;
    //     for (int i = 0; i < this->Nx; i++)
    //     {
    //         this->U_plus_three_quarter[time][i] = this->U_plus_current[i];
    //         this->U_minus_three_quarter[time][i] = this->U_minus_current[i];
    //     }
    // }

    // show the percent of completed computation

    for (int percent = 1; percent <= 10; ++percent)
    {
        if (n < floor(this->Nt * percent / 10))
        {
            break;
        }
        else if (n == floor(this->Nt * percent / 10))
        {
            cout << "The computation is completed " << 10 * percent << "%!" << endl;
            break;
        }
        else if (n == this->Nt - 2)
        {
            cout << "The computation is completed 100%!" << endl;
            break;
        }
    }
}

// void Solver::Fractional_step_method_scheme(int order_RK)
// {
//     this->copy_value_current_to_save();    // save the value of U_current to U_save
//     this->Upwind_scheme_non_source_term(); // compute U_next by Upwind scheme with non-source
//     this->copy_value_next_to_current();    // copy the value of U_next to U_current to prepare for ODE solver
//     this->Runge_Kutta_method(order_RK);    // compute U_next by RK method from U_current
//     this->copy_value_save_to_current();    // return the exact value for U_current at the current step
// }

void Solver::Fractional_step_method_limiter_Runge_Kutta_scheme()
{
    this->copy_value_current_to_save();            // save the value of U_current to U_save
    this->High_resoluion_method_non_source_term(); // compute U_next by high resoluion method using MC limiter with non-source
    this->copy_value_next_to_current();            // copy the value of U_next to U_current to prepare for ODE solver
    this->Runge_Kutta_method(2);                   // compute U_next by RK2 method from U_current
    this->copy_value_save_to_current();            // return the exact value for U_current at the current step
}

void Solver::MacCormack_scheme()
{
    // Save the value of density at current step before computing
    this->copy_value_current_to_save();

    this->Upwind_scheme();
    this->copy_value_next_to_current();
    this->Downwind_scheme();
    for (int i = 0; i < this->Nx; i++)
    {
        this->U_plus_next[i] = 0.5 * (this->U_plus_current_save[i] + this->U_plus_next[i]);
        this->U_minus_next[i] = 0.5 * (this->U_minus_current_save[i] + this->U_minus_next[i]);
    }

    this->copy_value_save_to_current();
}

void Solver::Downwind_scheme()
{
    this->compute_source_term();

    // compute U_plus, U_minus
    double U_plus_ip, U_plus_i, U_minus_i, U_minus_im;

    for (int i = 0; i < this->Nx; i++)
    {
        U_plus_ip = get_U_plus_current(i + 1);
        U_plus_i = get_U_plus_current(i);
        U_minus_i = get_U_minus_current(i);
        U_minus_im = get_U_minus_current(i - 1);

        this->U_plus_next[i] = U_plus_i - dt / dx * gamma * (U_plus_ip - U_plus_i) + dt * get_Source_plus(i);

        this->U_minus_next[i] = U_minus_i - dt / dx * (-gamma) * (U_minus_i - U_minus_im) + dt * get_Source_minus(i);
    }
}

void Solver::Upwind_scheme()
{
    this->compute_source_term();

    // compute U_plus, U_minus
    double U_plus_i, U_plus_im, U_minus_ip, U_minus_i;

    for (int i = 0; i < this->Nx; i++)
    {
        U_plus_i = get_U_plus_current(i);
        U_plus_im = get_U_plus_current(i - 1);
        U_minus_i = get_U_minus_current(i);
        U_minus_ip = get_U_minus_current(i + 1);

        this->U_plus_next[i] = U_plus_i - dt / dx * gamma * (U_plus_i - U_plus_im) + dt * get_Source_plus(i);

        this->U_minus_next[i] = U_minus_i - dt / dx * (-gamma) * (U_minus_ip - U_minus_i) + dt * get_Source_minus(i);
    }
}

// void Solver::Upwind_scheme_non_source_term()
// {
//     double slope_MC = 0;

//     // compute U_plus
//     for (int i = 1; i < this->Nx; i++)
//     {
//         this->U_plus_next[i] = U_plus_current[i] - dt / dx * gamma * (U_plus_current[i] - U_plus_current[i - 1]);
//     }
//     if (this->boundary_condition == PBC)
//     {
//         U_plus_next[0] = U_plus_next[this->Nx - 1];
//     }
//     else
//     {
//         exit(EXIT_FAILURE);
//     }

//     // compute U_minus
//     for (int i = 0; i < this->Nx - 1; i++)
//     {
//         this->U_minus_next[i] = U_minus_current[i] - dt / dx * (-gamma) * (U_minus_current[i + 1] - U_minus_current[i]);
//     }
//     if (this->boundary_condition == PBC)
//     {
//         U_minus_next[this->Nx - 1] = U_minus_next[0];
//     }
//     else
//     {
//         exit(EXIT_FAILURE);
//     }
// }

void Solver::High_resoluion_method_non_source_term()
{
    // compute U_plus, U_minus
    double U_plus_i, U_plus_im, U_minus_ip, U_minus_i;

    for (int i = 0; i < this->Nx; i++)
    {
        U_plus_i = get_U_plus_current(i);
        U_plus_im = get_U_plus_current(i - 1);
        U_minus_i = get_U_minus_current(i);
        U_minus_ip = get_U_minus_current(i + 1);

        this->U_plus_next[i] = U_plus_i - dt / dx * gamma * (U_plus_i - U_plus_im) - 0.5 * gamma * dt / dx * (dx - gamma * dt) * (limiter_plus(i) - limiter_plus(i - 1));

        this->U_minus_next[i] = U_minus_i - dt / dx * (-gamma) * (U_minus_ip - U_minus_i) - 0.5 * gamma * dt / dx * (dx - gamma * dt) * (limiter_minus(i + 1) - limiter_minus(i));
    }
}

void Solver::Quasi_Steady_Algorithm_scheme()
{
    this->compute_source_term();

    // compute U_plus, U_minus
    double U_plus_i, U_plus_im, U_minus_ip, U_minus_i;

    for (int i = 0; i < this->Nx; i++)
    {
        U_plus_i = get_U_plus_current(i);
        U_plus_im = get_U_plus_current(i - 1);
        U_minus_i = get_U_minus_current(i);
        U_minus_ip = get_U_minus_current(i + 1);

        this->U_plus_next[i] = U_plus_i - CFL * (U_plus_i - U_plus_im) + 0.5 * dt * (get_Source_plus(i) + get_Source_plus(i - 1));

        this->U_minus_next[i] = U_minus_i + CFL * (U_minus_ip - U_minus_i) + 0.5 * dt * (get_Source_minus(i + 1) + get_Source_minus(i));
    }
}

void Solver::Quasi_Steady_Algorithm_limiter_scheme()
{
    this->compute_source_term();

    // compute U_plus, U_minus
    double U_plus_i, U_plus_im, U_minus_ip, U_minus_i;

    for (int i = 0; i < this->Nx; i++)
    {
        U_plus_i = get_U_plus_current(i);
        U_plus_im = get_U_plus_current(i - 1);
        U_minus_i = get_U_minus_current(i);
        U_minus_ip = get_U_minus_current(i + 1);

        this->U_plus_next[i] = U_plus_i - dt / dx * gamma * (U_plus_i - U_plus_im) - 0.5 * gamma * dt / dx * (dx - gamma * dt) * (limiter_plus_quasi_steady(i) - limiter_plus_quasi_steady(i - 1)) + 0.5 * dt * (get_Source_plus(i) + get_Source_plus(i - 1));

        this->U_minus_next[i] = U_minus_i - dt / dx * (-gamma) * (U_minus_ip - U_minus_i) - 0.5 * gamma * dt / dx * (dx - gamma * dt) * (limiter_minus_quasi_steady(i + 1) - limiter_minus_quasi_steady(i)) + 0.5 * dt * (get_Source_minus(i + 1) + get_Source_minus(i));
    }
}

void Solver::Quasi_Steady_Algorithm_Center_scheme()
{
    this->compute_source_term();

    // compute U_plus, U_minus
    double U_plus_ip, U_plus_i, U_plus_im, U_plus_im2, Source_plus_ip, Source_plus_i, Source_plus_im, Source_plus_im2;
    double U_minus_im, U_minus_i, U_minus_ip, U_minus_ip2, Source_minus_im, Source_minus_i, Source_minus_ip, Source_minus_ip2;

    for (int i = 0; i < this->Nx; i++)
    {
        U_plus_ip = get_U_plus_current(i + 1);
        U_plus_i = get_U_plus_current(i);
        U_plus_im = get_U_plus_current(i - 1);
        U_plus_im2 = get_U_plus_current(i - 2);

        Source_plus_ip = get_Source_plus(i + 1);
        Source_plus_i = get_Source_plus(i);
        Source_plus_im = get_Source_plus(i - 1);
        Source_plus_im2 = get_Source_plus(i - 2);

        U_minus_im = get_U_minus_current(i - 1);
        U_minus_i = get_U_minus_current(i);
        U_minus_ip = get_U_minus_current(i + 1);
        U_minus_ip2 = get_U_minus_current(i + 2);

        Source_minus_im = get_Source_minus(i - 1);
        Source_minus_i = get_Source_minus(i);
        Source_minus_ip = get_Source_minus(i + 1);
        Source_minus_ip2 = get_Source_minus(i + 2);

        this->U_plus_next[i] = U_plus_i - 0.25 * CFL * (U_plus_ip + 3 * U_plus_i - 5 * U_plus_im + U_plus_im2) + 0.25 * CFL * CFL * (U_plus_ip - U_plus_i - U_plus_im + U_plus_im2) + 0.125 * dt * (Source_plus_ip + 5 * Source_plus_i + 3 * Source_plus_im - Source_plus_im2) - 0.125 * CFL * dt * (Source_plus_ip + Source_plus_i - Source_plus_im - Source_plus_im2);

        this->U_minus_next[i] = U_minus_i - 0.25 * CFL * (U_minus_im + 3 * U_minus_i - 5 * U_minus_ip + U_minus_ip2) + 0.25 * CFL * CFL * (U_minus_im - U_minus_i - U_minus_ip + U_minus_ip2) + 0.125 * dt * (Source_minus_im + 5 * Source_minus_i + 3 * Source_minus_ip - Source_minus_ip2) - 0.125 * CFL * dt * (Source_minus_im + Source_minus_i - Source_minus_ip - Source_minus_ip2) ;
    }
}

void Solver::Quasi_Steady_Algorithm_Beam_Warming_scheme()
{
    this->compute_source_term();

    // compute U_plus, U_minus
    double U_plus_i, U_plus_im, U_plus_im2, Source_plus_i, Source_plus_im, Source_plus_im2;
    double U_minus_i, U_minus_ip, U_minus_ip2, Source_minus_i, Source_minus_ip, Source_minus_ip2;

    for (int i = 0; i < this->Nx; i++)
    {
        U_plus_i = get_U_plus_current(i);
        U_plus_im = get_U_plus_current(i - 1);
        U_plus_im2 = get_U_plus_current(i - 2);

        Source_plus_i = get_Source_plus(i);
        Source_plus_im = get_Source_plus(i - 1);
        Source_plus_im2 = get_Source_plus(i - 2);

        U_minus_i = get_U_minus_current(i);
        U_minus_ip = get_U_minus_current(i + 1);
        U_minus_ip2 = get_U_minus_current(i + 2);

        Source_minus_i = get_Source_minus(i);
        Source_minus_ip = get_Source_minus(i + 1);
        Source_minus_ip2 = get_Source_minus(i + 2);

        this->U_plus_next[i] = U_plus_i - 0.5 * CFL * (3 * U_plus_i - 4 * U_plus_im + U_plus_im2) + 0.5 * CFL * CFL * (U_plus_i - 2 * U_plus_im + U_plus_im2) + 0.25 * dt * (3 * Source_plus_i + 2 * Source_plus_im - Source_plus_im2) - 0.25 * CFL * dt * (Source_plus_i - Source_plus_im2) ;

        this->U_minus_next[i] = U_minus_i - 0.5 * CFL * (3 * U_minus_i - 4 * U_minus_ip + U_minus_ip2) + 0.5 * CFL * CFL * (U_minus_i - 2 * U_minus_ip + U_minus_ip2) + 0.25 * dt * (3 * Source_minus_i + 2 * Source_minus_ip - Source_minus_ip2) - 0.25 * CFL * dt * (Source_minus_i - Source_minus_ip2) ;
    }
}

void Solver::Quasi_Steady_Algorithm_Lax_Wendroof_scheme()
{
    this->compute_source_term();

    // compute U_plus, U_minus
    double U_plus_i, U_plus_im, U_plus_ip, Source_plus_i, Source_plus_im, Source_plus_ip;
    double U_minus_i, U_minus_im, U_minus_ip, Source_minus_i, Source_minus_im, Source_minus_ip;

    for (int i = 0; i < this->Nx; i++)
    {
        U_plus_i = get_U_plus_current(i);
        U_plus_im = get_U_plus_current(i - 1);
        U_plus_ip = get_U_plus_current(i + 1);

        Source_plus_i = get_Source_plus(i);
        Source_plus_im = get_Source_plus(i - 1);
        Source_plus_ip = get_Source_plus(i + 1);

        U_minus_i = get_U_minus_current(i);
        U_minus_ip = get_U_minus_current(i + 1);
        U_minus_im = get_U_minus_current(i - 1);

        Source_minus_i = get_Source_minus(i);
        Source_minus_ip = get_Source_minus(i + 1);
        Source_minus_im = get_Source_minus(i - 1);

        this->U_plus_next[i] = U_plus_i - 0.5 * CFL * (U_plus_ip - U_plus_im) + 0.5 * CFL * CFL * (U_plus_ip - 2 * U_plus_i + U_plus_im) + 0.25 * dt * (Source_plus_ip + 2 * Source_plus_i + Source_plus_im) - 0.25 * CFL * dt * (Source_plus_ip - Source_plus_im) ;

        this->U_minus_next[i] = U_minus_i - 0.5 * CFL * (U_minus_im - U_minus_ip) + 0.5 * CFL * CFL * (U_minus_im - 2 * U_minus_i + U_minus_ip) + 0.25 * dt * (Source_minus_im + 2 * Source_minus_i + Source_minus_ip) - 0.25 * CFL * dt * (Source_minus_im - Source_minus_ip) ;
    }
}

void Solver::F_Wave_Algorithm_scheme()
{
    this->copy_value_current_to_save();

    for (int i = 0; i < this->Nx - 1; i++)
    {
        U_plus_current[i] = 0.5 * (U_plus_current_save[i] + U_plus_current_save[i + 1]);
        U_minus_current[i] = 0.5 * (U_minus_current_save[i] + U_minus_current_save[i + 1]);
    }
    U_plus_current[this->Nx - 1] = U_plus_current[0];
    U_minus_current[this->Nx - 1] = U_minus_current[0];

    this->compute_source_term();

    // compute U_plus, U_minus
    double U_plus_ip, U_plus_i, U_plus_im, U_minus_ip, U_minus_i, U_minus_im, Source_plus_i, Source_plus_im, Source_minus_i, Source_minus_im;
    for (int i = 0; i < this->Nx; i++)
    {
        U_plus_ip = get_U_plus_current(i + 1);
        U_plus_i = get_U_plus_current(i);
        U_plus_im = get_U_plus_current(i - 1);
        U_minus_im = get_U_minus_current(i - 1);
        U_minus_i = get_U_minus_current(i);
        U_minus_ip = get_U_minus_current(i + 1);
        Source_plus_i = get_Source_plus(i);
        Source_plus_im = get_Source_plus(i - 1);
        Source_minus_i = get_Source_minus(i);
        Source_minus_im = get_Source_minus(i - 1);

        this->U_plus_next[i] = U_plus_i - 0.5 * gamma * dt / dx * (U_plus_ip - U_plus_im) + 0.5 * gamma * gamma * dt * dt / dx / dx * (U_plus_ip - 2 * U_plus_i + U_plus_im) + 0.5 * dt * (Source_plus_i + Source_plus_im) - 0.5 * gamma * dt * dt / dx * (Source_plus_i - Source_plus_im);

        this->U_minus_next[i] = U_minus_i - 0.5 * (-gamma) * dt / dx * (U_minus_ip - U_minus_im) + 0.5 * gamma * gamma * dt * dt / dx / dx * (U_minus_ip - 2 * U_minus_i + U_minus_im) + 0.5 * dt * (Source_minus_i + Source_minus_im) - 0.5 * (-gamma) * dt * dt / dx * (Source_minus_i - Source_minus_im);
    }

    this->copy_value_save_to_current();
}

double Solver::limiter_plus(int i)
{
    if (this->slope_limiter == Non_limiter)
    {
        return 0;
    }
    else
    {
        double U_ip, U_i, U_im, slope_center, slope_up, slope_down;

        U_ip = get_U_plus_current(i + 1);
        U_i = get_U_plus_current(i);
        U_im = get_U_plus_current(i - 1);

        slope_center = (U_ip - U_im) / (2 * this->dx);
        slope_up = (U_i - U_im) / (this->dx);
        slope_down = (U_ip - U_i) / (this->dx);

        return compute_slope(slope_center, slope_up, slope_down, "plus");
    }
}

double Solver::limiter_minus(int i)
{
    if (this->slope_limiter == Non_limiter)
    {
        return 0;
    }
    else
    {
        double U_ip, U_i, U_im, slope_center, slope_up, slope_down;

        U_ip = get_U_minus_current(i + 1);
        U_i = get_U_minus_current(i);
        U_im = get_U_minus_current(i - 1);

        slope_center = (U_ip - U_im) / (2 * this->dx);
        slope_up = (U_i - U_im) / (this->dx);
        slope_down = (U_ip - U_i) / (this->dx);

        return compute_slope(slope_center, slope_up, slope_down, "minus");
    }
}

double Solver::limiter_plus_quasi_steady(int i)
{
    double U_ip, U_i, U_im, delta_ip, delta_i, delta_im, slope_center, slope_up, slope_down;

    U_ip = get_U_plus_current(i + 1);
    delta_ip = 0.5 * dx / gamma * get_Source_plus(i + 1);
    U_i = get_U_plus_current(i);
    delta_i = 0.5 * dx / gamma * get_Source_plus(i);
    U_im = get_U_plus_current(i - 1);
    delta_im = 0.5 * dx / gamma * get_Source_plus(i - 1);

    slope_center = (U_ip - U_im - delta_ip - 2 * delta_i - delta_im) / (2 * this->dx);
    slope_up = (U_i - U_im - delta_i - delta_im) / (this->dx);
    slope_down = (U_ip - U_i - delta_ip - delta_i) / (this->dx);

    return compute_slope(slope_center, slope_up, slope_down, "plus");
}

double Solver::limiter_minus_quasi_steady(int i)
{
    double U_ip, U_i, U_im, delta_ip, delta_i, delta_im, slope_center, slope_up, slope_down;

    U_ip = get_U_minus_current(i + 1);
    delta_ip = -0.5 * dx / gamma * get_Source_minus(i + 1);
    U_i = get_U_minus_current(i);
    delta_i = -0.5 * dx / gamma * get_Source_minus(i);
    U_im = get_U_minus_current(i - 1);
    delta_im = -0.5 * dx / gamma * get_Source_minus(i - 1);

    slope_center = (U_ip - U_im - delta_ip - 2 * delta_i - delta_im) / (2 * this->dx);
    slope_up = (U_i - U_im - delta_i - delta_im) / (this->dx);
    slope_down = (U_ip - U_i - delta_ip - delta_i) / (this->dx);

    return compute_slope(slope_center, slope_up, slope_down, "minus");
}

double Solver::compute_slope(double slope_center, double slope_up, double slope_down, string plus_or_minus)
{
    switch (this->slope_limiter)
    {
    case Centered_slope:
        return slope_center;
        break;

    case Upwind_slope:
        return slope_up;

    case Downwind_slope:
        return slope_down;

    case Beam_Warming:
        if (plus_or_minus == "plus")
        {
            return slope_up;
        }
        else if (plus_or_minus == "minus")
        {
            return slope_down;
        }
        else
        {
            cerr << "There exits error in the limiter function" << endl;
            exit(EXIT_FAILURE);
        }

    case Lax_Wendroof:
        if (plus_or_minus == "plus")
        {
            return slope_down;
        }
        else if (plus_or_minus == "minus")
        {
            return slope_up;
        }
        else
        {
            cerr << "There exits error in the limiter function" << endl;
            exit(EXIT_FAILURE);
        }

    case Minmod:
        return minmod(slope_up, slope_down);

    case Superbee:
        return maxmod(minmod(slope_down, 2 * slope_up), minmod(2 * slope_down, slope_up));

    case MC_limiter:
        return minmod(slope_center, 2 * slope_up, 2 * slope_down);

    default:
        cerr << "There exits error in the limiter function" << endl;
        exit(EXIT_FAILURE);
    }
}

void Solver::check_tolerance_condition(int n)
{
    double L1_error_2 = 0;
    int time_n = int(n * this->dt);
    int time_n_plus = int((n + 1) * this->dt);

    if (time_n < time_n_plus)
    {
        L1_error_2 = this->compute_discrete_error();
        L1_error[time_n_plus - 1] = L1_error_2;
        if (satisfied_tolerance == false && L1_error_2 < this->tolerance && this->tolerance < 0.5)
        {
            satisfied_tolerance = true;

            if(this->initial_amplitude < 10.0)
            {
                this->Nt = int(min_of(this->Nt, 3.0 * n, n + 14000/this->dt));
            }
            else if(this->initial_amplitude < 20.0)
            {
                this->Nt = int(min_of(this->Nt, 2.5 * n, n + 10000/this->dt));
            }
            else if(this->initial_amplitude < 30.0)
            {
                this->Nt = int(min_of(this->Nt, 2.0 * n, n + 8000/this->dt));
            }
            else
            {
                this->Nt = int(min_of(this->Nt, 1.8 * n, n + 4000/this->dt));
            }

            this->final_time = ceil(this->Nt * this->dt);
        }

        //check again the tolerance condition for the case trainsient solution E(t) < 10^-14
        if(satisfied_tolerance == true && L1_error_2 > this->tolerance)
        {
            satisfied_tolerance = false;
            this->Nt = this->Nt_input;
            this->final_time = this->final_time_input;
        }

        //check non-convergence
        if(non_convergence == true && L1_error_2 < non_convergence_error_min)
        {
            this->non_convergence = false;
        }
        if(non_convergence == true && time_n_plus > non_convergence_stop_time && this->tolerance < 0.5) //stop the simulation
        {
            this->non_convergence = false;
            this->Nt = n + 2*number_of_final_step; 
            this->final_time = ceil(this->Nt * this->dt);
        }
    }
}

void Solver::compute_source_term()
{
    if (parallel == false)
    {
        this->compute_turning_rate();
    }
    else
    {
        this->compute_turning_rate_parallel();
    }

    for (int i = 0; i < this->Nx; i++)
    {
        this->Source_term_plus[i] = -Lambda_plus[i] * U_plus_current[i] + Lambda_minus[i] * U_minus_current[i];
        this->Source_term_minus[i] = -this->Source_term_plus[i];
    }
}

void Solver::copy_value_next_to_current()
{

    this->U_plus_current = this->U_plus_next;
    this->U_minus_current = this->U_minus_next;
}

void Solver::copy_value_save_to_current()
{

    this->U_plus_current = this->U_plus_current_save;
    this->U_minus_current = this->U_minus_current_save;
}

void Solver::copy_value_current_to_save()
{

    this->U_plus_current_save = this->U_plus_current;
    this->U_minus_current_save = this->U_minus_current;
}

void Solver::Runge_Kutta_method(int order)
{
    if (order == 2)
    {
        this->U_plus_temp = this->U_plus_current;
        this->U_minus_temp = this->U_minus_current;

        this->compute_source_term();

        for (int i = 0; i < this->Nx; i++)
        {
            this->U_plus_current[i] = this->U_plus_temp[i] + 0.5 * this->dt * this->Source_term_plus[i];
            this->U_minus_current[i] = this->U_minus_temp[i] + 0.5 * this->dt * this->Source_term_minus[i];
        }

        this->compute_source_term();

        // compute next step
        for (int i = 0; i < this->Nx; i++)
        {
            this->U_plus_next[i] = this->U_plus_temp[i] + this->dt * this->Source_term_plus[i];
            this->U_minus_next[i] = this->U_minus_temp[i] + this->dt * this->Source_term_minus[i];
        }
    }
    else if (order == 4)
    {
        this->U_plus_temp = this->U_plus_current;
        this->U_minus_temp = this->U_minus_current;

        // compute RK1
        this->compute_source_term();

        for (int i = 0; i < this->Nx; i++)
        {
            Runge_Kutta_1p[i] = this->dt * this->Source_term_plus[i];
            Runge_Kutta_1m[i] = this->dt * this->Source_term_minus[i];
        }

        // compute RK2
        for (int i = 0; i < this->Nx; i++)
        {
            this->U_plus_current[i] = this->U_plus_temp[i] + 0.5 * Runge_Kutta_1p[i];
            this->U_minus_current[i] = this->U_minus_temp[i] + 0.5 * Runge_Kutta_1m[i];
        }

        this->compute_source_term();

        for (int i = 0; i < this->Nx; i++)
        {
            Runge_Kutta_2p[i] = this->dt * this->Source_term_plus[i];
            Runge_Kutta_2m[i] = this->dt * this->Source_term_minus[i];
        }

        // compute RK3
        for (int i = 0; i < this->Nx; i++)
        {
            this->U_plus_current[i] = this->U_plus_temp[i] + 0.5 * Runge_Kutta_2p[i];
            this->U_minus_current[i] = this->U_minus_temp[i] + 0.5 * Runge_Kutta_2m[i];
        }

        this->compute_source_term();

        for (int i = 0; i < this->Nx; i++)
        {
            Runge_Kutta_3p[i] = this->dt * this->Source_term_plus[i];
            Runge_Kutta_3m[i] = this->dt * this->Source_term_minus[i];
        }

        // compute RK4
        for (int i = 0; i < this->Nx; i++)
        {
            this->U_plus_current[i] = this->U_plus_temp[i] + Runge_Kutta_3p[i];
            this->U_minus_current[i] = this->U_minus_temp[i] + Runge_Kutta_3m[i];
        }

        this->compute_source_term();

        for (int i = 0; i < this->Nx; i++)
        {
            Runge_Kutta_4p[i] = this->dt * this->Source_term_plus[i];
            Runge_Kutta_4m[i] = this->dt * this->Source_term_minus[i];
        }

        // compute next step
        for (int i = 0; i < this->Nx; i++)
        {
            this->U_plus_next[i] = this->U_plus_temp[i] + 1 / 6.0 * (Runge_Kutta_1p[i] + 2 * Runge_Kutta_2p[i] + 2 * Runge_Kutta_3p[i] + Runge_Kutta_4p[i]);
            this->U_minus_next[i] = this->U_minus_temp[i] + 1 / 6.0 * (Runge_Kutta_1m[i] + 2 * Runge_Kutta_2m[i] + 2 * Runge_Kutta_3m[i] + Runge_Kutta_4m[i]);
        }
    }
    else
    {
        cerr << "There is only Runge Kutta method order two or four!" << endl;
        exit(EXIT_FAILURE);
    }
}

void Solver::compute_turning_rate()
{
    double Y_r_p = 0; // y_r^+
    // double Y_r_m = 0;  // y_r^-
    double Y_a_p = 0; // y_a^+
    // double Y_a_m = 0;  // y_a^-
    double Y_al_p = 0; // y_al^+
    // double Y_al_m = 0; // y_al^-

    for (int i = 0; i < this->Nx; i++)
    {
        if (this->q_r != 0)
        {
            Y_r_p = this->q_r * compute_integral(Repulsion, i); // y_r^+
        }

        if (this->q_a != 0)
        {
            Y_a_p = this->q_a * compute_integral(Attraction, i); // y_a^+
        }

        if (this->q_al != 0)
        {
            Y_al_p = this->q_al * compute_integral(Alignment, i); // y_al^+
        }

        this->Lambda_plus[i] = lambda_1 + lambda_2 * (0.5 + 0.5 * tanh(+Y_r_p - Y_a_p + Y_al_p - y_0));
        this->Lambda_minus[i] = lambda_1 + lambda_2 * (0.5 + 0.5 * tanh(-Y_r_p + Y_a_p - Y_al_p - y_0));
    }
}

void Solver::compute_turning_rate_parallel()
{
    if (this->number_of_threads > 0 && this->number_of_threads != 4)
    {
        omp_set_num_threads(this->number_of_threads);
        omp_set_nested(true);
#pragma omp parallel for
        for (int i = 0; i < this->Nx; i++)
        {
            double Y_r_p = 0.0;
            double Y_a_p = 0.0;
            double Y_al_p = 0.0;

            if (this->q_r != 0)
            {
                Y_r_p = this->q_r * compute_integral(Repulsion, i); // y_r^+
            }

            if (this->q_a != 0)
            {
                Y_a_p = this->q_a * compute_integral(Attraction, i); // y_a^+
            }

            if (this->q_al != 0)
            {
                Y_al_p = this->q_al * compute_integral(Alignment, i); // y_al^+
            }

            this->Lambda_plus[i] = lambda_1 + lambda_2 * (0.5 + 0.5 * tanh(+Y_r_p - Y_a_p + Y_al_p - y_0));
            this->Lambda_minus[i] = lambda_1 + lambda_2 * (0.5 + 0.5 * tanh(-Y_r_p + Y_a_p - Y_al_p - y_0));
        }
    }
    else
    {
        // use 4 parallel sections
        int parallel_section_1 = int(this->Nx * 1 / 4);
        int parallel_section_2 = int(this->Nx * 2 / 4);
        int parallel_section_3 = int(this->Nx * 3 / 4);
        int parallel_section_4 = this->Nx;

        omp_set_num_threads(4);
#pragma omp parallel sections
        {
#pragma omp section // parallel section 1
            {
                for (int i = 0; i < parallel_section_1; ++i)
                {
                    double Y_r_p = 0.0;
                    double Y_a_p = 0.0;
                    double Y_al_p = 0.0;

                    if (this->q_r != 0)
                    {
                        Y_r_p = this->q_r * compute_integral(Repulsion, i); // y_r^+
                    }

                    if (this->q_a != 0)
                    {
                        Y_a_p = this->q_a * compute_integral(Attraction, i); // y_a^+
                    }

                    if (this->q_al != 0)
                    {
                        Y_al_p = this->q_al * compute_integral(Alignment, i); // y_al^+
                    }

                    this->Lambda_plus[i] = lambda_1 + lambda_2 * (0.5 + 0.5 * tanh(+Y_r_p - Y_a_p + Y_al_p - y_0));
                    this->Lambda_minus[i] = lambda_1 + lambda_2 * (0.5 + 0.5 * tanh(-Y_r_p + Y_a_p - Y_al_p - y_0));
                }
            }
#pragma omp section // parallel section 2
            {
                for (int i = parallel_section_1; i < parallel_section_2; ++i)
                {
                    double Y_r_p = 0.0;
                    double Y_a_p = 0.0;
                    double Y_al_p = 0.0;

                    if (this->q_r != 0)
                    {
                        Y_r_p = this->q_r * compute_integral(Repulsion, i); // y_r^+
                    }

                    if (this->q_a != 0)
                    {
                        Y_a_p = this->q_a * compute_integral(Attraction, i); // y_a^+
                    }

                    if (this->q_al != 0)
                    {
                        Y_al_p = this->q_al * compute_integral(Alignment, i); // y_al^+
                    }

                    this->Lambda_plus[i] = lambda_1 + lambda_2 * (0.5 + 0.5 * tanh(+Y_r_p - Y_a_p + Y_al_p - y_0));
                    this->Lambda_minus[i] = lambda_1 + lambda_2 * (0.5 + 0.5 * tanh(-Y_r_p + Y_a_p - Y_al_p - y_0));
                }
            }

#pragma omp section // parallel section 3
            {
                for (int i = parallel_section_2; i < parallel_section_3; ++i)
                {
                    double Y_r_p = 0.0;
                    double Y_a_p = 0.0;
                    double Y_al_p = 0.0;

                    if (this->q_r != 0)
                    {
                        Y_r_p = this->q_r * compute_integral(Repulsion, i); // y_r^+
                    }

                    if (this->q_a != 0)
                    {
                        Y_a_p = this->q_a * compute_integral(Attraction, i); // y_a^+
                    }

                    if (this->q_al != 0)
                    {
                        Y_al_p = this->q_al * compute_integral(Alignment, i); // y_al^+
                    }

                    this->Lambda_plus[i] = lambda_1 + lambda_2 * (0.5 + 0.5 * tanh(+Y_r_p - Y_a_p + Y_al_p - y_0));
                    this->Lambda_minus[i] = lambda_1 + lambda_2 * (0.5 + 0.5 * tanh(-Y_r_p + Y_a_p - Y_al_p - y_0));
                }
            }

#pragma omp section // parallel section 4
            {
                for (int i = parallel_section_3; i < parallel_section_4; ++i)
                {
                    double Y_r_p = 0.0;
                    double Y_a_p = 0.0;
                    double Y_al_p = 0.0;

                    if (this->q_r != 0)
                    {
                        Y_r_p = this->q_r * compute_integral(Repulsion, i); // y_r^+
                    }

                    if (this->q_a != 0)
                    {
                        Y_a_p = this->q_a * compute_integral(Attraction, i); // y_a^+
                    }

                    if (this->q_al != 0)
                    {
                        Y_al_p = this->q_al * compute_integral(Alignment, i); // y_al^+
                    }

                    this->Lambda_plus[i] = lambda_1 + lambda_2 * (0.5 + 0.5 * tanh(+Y_r_p - Y_a_p + Y_al_p - y_0));
                    this->Lambda_minus[i] = lambda_1 + lambda_2 * (0.5 + 0.5 * tanh(-Y_r_p + Y_a_p - Y_al_p - y_0));
                }
            }
        }
    }
}

void Solver::check_all_data()
{
    Computation_data::check_computation_data();

    if (this->scheme == Upwind || this->scheme == MacCormack)
    {
        // check the convective instability (the CFL condition)
        double CFL_value = this->gamma * this->dt / this->dx;
        if (CFL_value > 1.0)
        {
            cerr << "CFL value: " << CFL_value << endl;
            cerr << "The CFL condition does not satisfy. Please make sure that CFL condition satisfies!" << endl;
            exit(EXIT_FAILURE);
        }

        // check the relaxation instability
        if (this->dt > 0.038)
        {
            cerr << "The time step: " << this->dt << endl;
            cerr << "The time step is not enough to ensure that there is no relaxation instability. Please check again!" << endl;
            exit(EXIT_FAILURE);
        }
    }

    // cout << "All data is OK!" << endl;
}

double Solver::compute_integral(Interaction interaction, int i)
{
    int Ns = 0; // The number of the step

    if (interaction == Repulsion)
    {
        Ns = Nr;
    }
    else if (interaction == Attraction)
    {
        Ns = Na;
    }
    else if (interaction == Alignment)
    {
        Ns = Nal;
    }

    double result = 0;
    double result_0 = compute_interaction_kernel_extension(interaction, i, 0);
    double result_end = compute_interaction_kernel_extension(interaction, i, Ns - 1);
    double result_odd = 0;
    double result_even = 0;
    bool flag = true;

    for (int s = 1; s <= Ns - 2; s++)
    {
        if (flag == true)
        {
            result_odd += compute_interaction_kernel_extension(interaction, i, s);
        }
        else
        {
            result_even += compute_interaction_kernel_extension(interaction, i, s);
        }

        flag = !flag;
    }

    result = this->dx / 3 * (result_0 + 4 * result_odd + 2 * result_even + result_end);
    return result;
}

double Solver::compute_interaction_kernel_extension(Interaction interaction, int i, int s)
{
    double s_j = 0; // the spatial region for attractive, alignment or repulsive
    double m_j = 0; // Width of attraction kernel, alignment kernel or repulsion kernel

    if (interaction == Repulsion)
    {
        s_j = this->s_r;
        m_j = this->m_r;
    }
    else if (interaction == Attraction)
    {
        s_j = this->s_a;
        m_j = this->m_a;
    }
    else if (interaction == Alignment)
    {
        s_j = this->s_al;
        m_j = this->m_al;
    }

    double K_j = 1 / sqrt(2 * M_PI * m_j * m_j) * exp(-(s * this->dx - s_j) * (s * this->dx - s_j) / (2 * m_j * m_j));

    double result = 0;
    if (this->model == M1)
    {
        if (interaction == Repulsion || interaction == Attraction)
        {
            result = get_U_plus_current(i + s) - get_U_plus_current(i - s) + get_U_minus_current(i + s) - get_U_minus_current(i - s);
        }
        else if (interaction == Alignment)
        {
            result = get_U_minus_current(i + s) - get_U_plus_current(i - s);
        }
    }
    else
    {
        cerr << "Please code for another model!" << endl;
        exit(EXIT_FAILURE);
    }
    return K_j * result;
}

double Solver::get_U_plus_current(int position)
{
    /*
    if ((position < -Nx) || (position > 2 * Nx - 1))
    {
        cerr << "The position is not correct!" << endl;
        exit(EXIT_FAILURE);
    }

    if ((n < 0) || (n > Nt - 1))
    {
        cerr << "The time is not correct!" << endl;
        exit(EXIT_FAILURE);
    }
    */

    if (boundary_condition == PBC)
    {
        if ((position >= 0) && (position < (Nx - 1)))
        {
            return this->U_plus_current[position];
        }

        if (position >= (Nx - 1))
        {
            return this->U_plus_current[position - (Nx - 1)];
        }

        if (position < 0)
        {
            return this->U_plus_current[(Nx - 1) + position];
        }
    }
    else
    {
        cerr << "Please code for another boundary condition!" << endl;
        exit(EXIT_FAILURE);
    }

    return 0;
}

double Solver::get_U_minus_current(int position)
{
    /*
    if ((position < -Nx) || (position > 2 * Nx - 1))
    {
        cerr << "The position is not correct!" << endl;
        exit(EXIT_FAILURE);
    }

    if ((n < 0) || (n > Nt - 1))
    {
        cerr << "The time is not correct!" << endl;
        exit(EXIT_FAILURE);
    }
    */

    if (boundary_condition == PBC)
    {
        if ((position >= 0) && (position < (Nx - 1)))
        {
            return this->U_minus_current[position];
        }

        if (position >= (Nx - 1))
        {
            return this->U_minus_current[position - (Nx - 1)];
        }

        if (position < 0)
        {
            return this->U_minus_current[(Nx - 1) + position];
        }
    }
    else
    {
        cerr << "Please code for another boundary condition!" << endl;
        exit(EXIT_FAILURE);
    }

    return 0;
}

double Solver::get_Lambda_plus(int position)
{
    /*
    if ((position < -Nx) || (position > 2 * Nx - 1))
    {
        cerr << "The position is not correct!" << endl;
        exit(EXIT_FAILURE);
    }
    */

    if (boundary_condition == PBC)
    {
        if ((position >= 0) && (position < (Nx - 1)))
        {
            return this->Lambda_plus[position];
        }

        if (position >= (Nx - 1))
        {
            return this->Lambda_plus[position - (Nx - 1)];
        }

        if (position < 0)
        {
            return this->Lambda_plus[(Nx - 1) + position];
        }
    }
    else
    {
        cerr << "Please code for another boundary condition!" << endl;
        exit(EXIT_FAILURE);
    }

    return 0;
}

double Solver::get_Lambda_minus(int position)
{
    if (boundary_condition == PBC)
    {
        if ((position >= 0) && (position < (Nx - 1)))
        {
            return this->Lambda_minus[position];
        }

        if (position >= (Nx - 1))
        {
            return this->Lambda_minus[position - (Nx - 1)];
        }

        if (position < 0)
        {
            return this->Lambda_minus[(Nx - 1) + position];
        }
    }
    else
    {
        cerr << "Please code for another boundary condition!" << endl;
        exit(EXIT_FAILURE);
    }

    return 0;
}

double Solver::get_Source_minus(int position)
{
    if (boundary_condition == PBC)
    {
        if ((position >= 0) && (position < (Nx - 1)))
        {
            return this->Source_term_minus[position];
        }

        if (position >= (Nx - 1))
        {
            return this->Source_term_minus[position - (Nx - 1)];
        }

        if (position < 0)
        {
            return this->Source_term_minus[(Nx - 1) + position];
        }
    }
    else
    {
        cerr << "Please code for another boundary condition!" << endl;
        exit(EXIT_FAILURE);
    }

    return 0;
}

double Solver::get_Source_plus(int position)
{
    if (boundary_condition == PBC)
    {
        if ((position >= 0) && (position < (Nx - 1)))
        {
            return this->Source_term_plus[position];
        }

        if (position >= (Nx - 1))
        {
            return this->Source_term_plus[position - (Nx - 1)];
        }

        if (position < 0)
        {
            return this->Source_term_plus[(Nx - 1) + position];
        }
    }
    else
    {
        cerr << "Please code for another boundary condition!" << endl;
        exit(EXIT_FAILURE);
    }

    return 0;
}

void Solver::create_objects()
{
    //this->three_quarter = int((this->Nt - 1) * 3 / 4);
    this->CFL = this->gamma * this->dt / this->dx;

    this->Nt_input = this->Nt;
    this->final_time_input = this->final_time;

    Resize(this->U_plus_stable, this->number_of_final_step, this->Nx);
    Resize(this->U_minus_stable, this->number_of_final_step, this->Nx);

    //Resize(this->U_plus_three_quarter, this->number_of_final_step, this->Nx);
    //Resize(this->U_minus_three_quarter, this->number_of_final_step, this->Nx);

    Resize(this->U_plus_current, this->Nx);
    Resize(this->U_minus_current, this->Nx);
    Resize(this->U_plus_next, this->Nx);
    Resize(this->U_minus_next, this->Nx);

    // Resize(this->Total_density_current, this->Nx);

    Resize(this->Lambda_plus, this->Nx);
    Resize(this->Lambda_minus, this->Nx);

    Resize(this->U_plus_temp, this->Nx);
    Resize(this->U_minus_temp, this->Nx);

    Resize(this->U_plus_temp, this->Nx);
    Resize(this->U_minus_temp, this->Nx);

    Resize(this->U_plus_current_save, this->Nx);
    Resize(this->U_minus_current_save, this->Nx);

    Resize(this->Source_term_plus, this->Nx);
    Resize(this->Source_term_minus, this->Nx);

    Resize(this->Runge_Kutta_1p, this->Nx);
    Resize(this->Runge_Kutta_2p, this->Nx);
    Resize(this->Runge_Kutta_3p, this->Nx);
    Resize(this->Runge_Kutta_4p, this->Nx);

    Resize(this->Runge_Kutta_1m, this->Nx);
    Resize(this->Runge_Kutta_2m, this->Nx);
    Resize(this->Runge_Kutta_3m, this->Nx);
    Resize(this->Runge_Kutta_4m, this->Nx);

    Resize(this->Total_desity_100, 100, this->Nx);
    Resize(this->L1_error, this->final_time);

    // cout << "Creating objects for solver is completed!" << endl;
}

void Solver::Solve_detail()
{
    double start_time = get_cpu_time_hour();

    cout << "2. Check the correctness of input data!" << endl;
    this->check_all_data();
    
    cout << "3. Delete all old output files!" << endl;
    this->delete_old_output();

    cout << "4. Create objects for solver before computations!" << endl;
    this->create_objects();

    cout << "5. Write the input data before computation!" << endl;
    this->write_input_computation_data();
    cout << "5. Solver by using " << Scheme_string[this->scheme] << " and " << Slope_string[this->slope_limiter] << endl;
    cout << "5. Magnitude of attraction q_a = " << this->q_a << "; The initial amplitude A = " << this->initial_amplitude << endl;
    cout << "5. The final time T =  " << this->final_time << "; The tolerance is " << this->tolerance << endl;
    cout << "5. The space step dx = " << this->dx << "; The time step dt = " << this->dt << endl;
    cout << "5. The output directory is " << this->output_directory << endl;

    cout << "6. Compute the solution: " << endl;
    this->compute_solution();

    cout << "7. Write the solution!" << endl;
    this->write_solution();
    
    cout << "8. Write the input data after computation!" << endl;
    double finish_time = get_cpu_time_hour();
    this->cpu_time_hour = finish_time - start_time;
    this->write_input_computation_data();
    cout << "8. Solver by using " << Scheme_string[this->scheme] << " and " << Slope_string[this->slope_limiter] << endl;
    cout << "8. Magnitude of attraction q_a = " << this->q_a << "; The initial amplitude A = " << this->initial_amplitude << endl;
    cout << "8. The final time T =  " << this->final_time << "; The tolerance is " << this->tolerance << endl;
    cout << "8. The space step dx = " << this->dx << "; The time step dt = " << this->dt << endl;
    cout << "8. The output directory is " << this->output_directory << endl;
    cout << "8. The CPU time is " << this->cpu_time_hour << " hours!" << endl;

    cout << "9. Finish successfully at " + get_local_time();
    cout << endl;
    cout << "----------------------------------------------------------------" << endl;
    cout << endl;
}

void Solver::Solve()
{
    this->parallel = false;
    this->number_of_threads = 1;
    this->Solve_detail();
}

void Solver::Solve_parallel()
{
    this->parallel = true;
    this->number_of_threads = 4; // default: use 4 CPUs
    this->Solve_detail();
}

void Solver::Solve_parallel(int number_threads)
{
    this->parallel = true;
    this->number_of_threads = number_threads;
    this->Solve_detail();
}

void Solver::write_solution()
{
    string filename_solution;

    // density of right-moving
    // filename_solution = this->output_directory + "/density_right_moving.txt";
    // Write_Data(this->U_plus, filename_solution.c_str(), time_scale, space_scale);

    // filename_solution = this->output_directory + "/density_right_moving_stability.txt";
    // Write_Data(this->U_plus_stable, filename_solution.c_str());

    //filename_solution = this->output_directory + "/density_right_moving_three_quarter_stable.txt";
    //Write_Data(this->U_plus_three_quarter, filename_solution.c_str());

    //filename_solution = this->output_directory + "/density_right_moving_three_quarter.txt";
    //Write_Data(this->U_plus_three_quarter.back(), filename_solution.c_str());

    filename_solution = this->output_directory + "/density_right_moving_finaltime.txt";
    Write_Data(this->U_plus_current, filename_solution.c_str());

    // density of left-moving
    // filename_solution = this->output_directory + "/density_left_moving.txt";
    // Write_Data(this->U_minus, filename_solution.c_str(), time_scale, space_scale);

    // filename_solution = this->output_directory + "/density_left_moving_stability.txt";
    // Write_Data(this->U_minus_stable, filename_solution.c_str());

    //filename_solution = this->output_directory + "/density_left_moving_three_quarter_stable.txt";
    //Write_Data(this->U_minus_three_quarter, filename_solution.c_str());

    //filename_solution = this->output_directory + "/density_left_moving_three_quarter.txt";
    //Write_Data(this->U_minus_three_quarter.back(), filename_solution.c_str());

    filename_solution = this->output_directory + "/density_left_moving_finaltime.txt";
    Write_Data(this->U_minus_current, filename_solution.c_str());

    // total density
    // filename_solution = this->output_directory + "/total_density.txt";
    // Write_Data(this->U_plus_stable, this->U_minus_stable, filename_solution.c_str(), time_scale, space_scale);

    filename_solution = this->output_directory + "/total_density_stability.txt";
    Write_Data(this->U_plus_stable, this->U_minus_stable, filename_solution.c_str());

    //filename_solution = this->output_directory + "/total_density_three_quarter_stable.txt";
    //Write_Data(this->U_plus_three_quarter, this->U_minus_three_quarter, filename_solution.c_str());

    //filename_solution = this->output_directory + "/total_density_three_quarter.txt";
    //Write_Data(this->U_plus_three_quarter.back(), this->U_minus_three_quarter.back(), filename_solution.c_str());

    filename_solution = this->output_directory + "/total_density_finaltime.txt";
    Write_Data(this->U_plus_current, this->U_minus_current, filename_solution.c_str());

    filename_solution = this->output_directory + "/total_density_100.txt";
    Write_Data(this->Total_desity_100, filename_solution.c_str());

    filename_solution = this->output_directory + "/L1_error.txt";
    Write_Data(this->L1_error, filename_solution.c_str());

    // cout << "Writting the solution is complete!" << endl;
}

void Solver::write_input_computation_data()
{
    // create the name of file
    string filename = this->output_directory + "/input_computation_data.txt";

    ofstream writefile(filename, ios::out | ios::app);

    // if the file can not be open then stop
    if (!writefile)
    {
        cerr << "Impossible to open the " << filename.c_str() << " file!" << endl;
        exit(EXIT_FAILURE);
    }

    writefile.precision(PRECISION);

    writefile << "### Writing time at " + get_local_time() << endl;

    writefile << "#####################################################################" << endl;
    writefile << "# Input data for computation" << endl;
    writefile << "#####################################################################" << endl;
    writefile << endl;

    writefile << "The left boundary of the domain <left_boundary>:: " << this->left_boundary << endl;
    writefile << endl;

    writefile << "The right boundary of the domain <right_boundary>:: " << this->right_boundary << endl;
    writefile << endl;

    writefile << "The step <dx>:: " << this->dx << endl;
    writefile << endl;

    writefile << "Final time <final_time>:: " << this->final_time << endl;
    writefile << endl;

    writefile << "Time step <dt>:: " << this->dt << endl;
    writefile << endl;

    writefile << "Tolerance of two adjacent steps <tolerance>:: " << this->tolerance << endl;
    writefile << endl;

    writefile << "The speed <gamma>:: " << this->gamma << endl;
    writefile << endl;

    writefile << "Base line turning rate <lambda_1>:: " << this->lambda_1 << endl;
    writefile << endl;

    writefile << "Bias turning rate <lambda_2>:: " << this->lambda_2 << endl;
    writefile << endl;

    writefile << "Magnitude of repulsion <q_r>:: " << this->q_r << endl;
    writefile << endl;

    writefile << "Magnitude of attraction <q_a>:: " << this->q_a << endl;
    writefile << endl;

    writefile << "Magnitude of alignment <q_al>:: " << this->q_al << endl;
    writefile << endl;

    writefile << "Shift of the turning function <y_0>:: " << this->y_0 << endl;
    writefile << endl;

    writefile << "The spatial region for attractive <s_a>:: " << this->s_a << endl;
    writefile << endl;

    writefile << "The spatial region for alignment <s_al>:: " << this->s_al << endl;
    writefile << endl;

    writefile << "The spatial region for repulsive <s_r>:: " << this->s_r << endl;
    writefile << endl;

    writefile << "Width of attraction kernel <m_a>:: " << this->m_a << endl;
    writefile << endl;

    writefile << "Width of alignment kernel <m_al>:: " << this->m_al << endl;
    writefile << endl;

    writefile << "Width of repulsion kernel <m_r>:: " << this->m_r << endl;
    writefile << endl;

    writefile << "The amplitude of the initial condition <initial_amplitude>:: " << this->initial_amplitude << endl;
    writefile << endl;

    writefile << "The model is used <model>:: " << Model_string[this->model] << endl;
    writefile << endl;

    writefile << "The scheme is used <scheme>:: " << Scheme_string[this->scheme] << endl;
    writefile << endl;

    writefile << "The slope limiter is used for FSM and QSA <slope_limiter>:: " << Slope_string[this->slope_limiter] << endl;
    writefile << endl;

    writefile << "Boundary condition <boundary_condition>:: " << Boundary_string[this->boundary_condition] << endl;
    writefile << endl;

    writefile << "The CPU time (hours):: " << this->cpu_time_hour << endl;
    writefile << endl;

    writefile << endl;
    writefile << endl;
    writefile.close();
}

VECTOR Solver::get_Total_density_final_time()
{
    VECTOR Total_density_final_time;
    Resize(Total_density_final_time, this->Nx);
    for (int i = 0; i < this->Nx; i++)
    {
        Total_density_final_time[i] = this->U_plus_current[i] + this->U_minus_current[i];
    }
    return Total_density_final_time;
}

/*
MATRIX Solver::get_U_plus()
{
    return this->U_plus;
}

MATRIX Solver::get_U_minus()
{
    return this->U_minus;
}
*/

int Solver::get_Nt()
{
    return this->Nt;
}

double Solver::get_dt()
{
    return this->dt;
}

int Solver::get_Nx()
{
    return this->Nx;
}

double Solver::get_dx()
{
    return this->dx;
}

int Solver::get_Nr()
{
    return this->Nr;
}

int Solver::get_Na()
{
    return this->Na;
}

int Solver::get_Nal()
{
    return this->Nal;
}

double Solver::get_CPU_time()
{
    return this->cpu_time_hour;
}

double Solver::get_initial_amplitude()
{
    return this->initial_amplitude;
}

double Solver::compute_discrete_error()
{
    double sum_error = 0;
    for (int i = 0; i < this->Nx; i++)
    {
        sum_error += abs(U_plus_next[i] + U_minus_next[i] - U_plus_current[i] - U_minus_current[i]);
    }
    return sum_error / this->Nx * (this->right_boundary - this->left_boundary);
}

void Solver::delete_old_output()
{
    const char *filename;
    string filename_old;

    filename_old = this->output_directory + "/input_computation_data.txt";
    filename = filename_old.c_str();
    remove(filename);

    filename_old = this->output_directory + "/density_right_moving_initial.txt";
    filename = filename_old.c_str();
    remove(filename);

    filename_old = this->output_directory + "/density_left_moving_initial.txt";
    filename = filename_old.c_str();
    remove(filename);

    filename_old = this->output_directory + "/total_density_initial.txt";
    filename = filename_old.c_str();
    remove(filename);

    filename_old = this->output_directory + "/density_right_moving_stability.txt";
    filename = filename_old.c_str();
    remove(filename);

    filename_old = this->output_directory + "/density_right_moving_three_quarter_stable.txt";
    filename = filename_old.c_str();
    remove(filename);

    filename_old = this->output_directory + "/density_right_moving_three_quarter.txt";
    filename = filename_old.c_str();
    remove(filename);

    filename_old = this->output_directory + "/density_right_moving_finaltime.txt";
    filename = filename_old.c_str();
    remove(filename);

    filename_old = this->output_directory + "/density_left_moving_stability.txt";
    filename = filename_old.c_str();
    remove(filename);

    filename_old = this->output_directory + "/density_left_moving_three_quarter_stable.txt";
    filename = filename_old.c_str();
    remove(filename);

    filename_old = this->output_directory + "/density_left_moving_three_quarter.txt";
    filename = filename_old.c_str();
    remove(filename);

    filename_old = this->output_directory + "/density_left_moving_finaltime.txt";
    filename = filename_old.c_str();
    remove(filename);

    filename_old = this->output_directory + "/total_density_stability.txt";
    filename = filename_old.c_str();
    remove(filename);

    filename_old = this->output_directory + "/total_density_three_quarter_stable.txt";
    filename = filename_old.c_str();
    remove(filename);

    filename_old = this->output_directory + "/total_density_three_quarter.txt";
    filename = filename_old.c_str();
    remove(filename);

    filename_old = this->output_directory + "/total_density_finaltime.txt";
    filename = filename_old.c_str();
    remove(filename);

    filename_old = this->output_directory + "/total_density_100.txt";
    filename = filename_old.c_str();
    remove(filename);

    filename_old = this->output_directory + "/L1_error.txt";
    filename = filename_old.c_str();
    remove(filename);

    // cout << "Deleted all old output files!" << endl;
}

void Solver::set_initial_amplitude(double new_initial_amplitude)
{
    this->initial_amplitude = new_initial_amplitude;
}

void Solver::set_scheme(Scheme new_scheme)
{
    this->scheme = new_scheme;
}

void Solver::set_slope_limiter(SlopeLimiter new_slope)
{
    this->slope_limiter = new_slope;
}

void Solver::set_q_a(double new_q_a)
{
    this->q_a = new_q_a;
}

void Solver::set_output_directory(string new_output_directory)
{
    this->output_directory = new_output_directory;
}

void Solver::set_final_time(double new_final_time)
{
    this->final_time = new_final_time;
    this->Nt = Get_Length(0, this->final_time, this->dt);
}

void Solver::set_dx(double new_dx)
{
    this->dx = new_dx;
    this->Nx = Get_Length(this->left_boundary, this->right_boundary, this->dx);
    this->Na = Get_Length(0, 2 * s_a, dx);
    this->Nr = Get_Length(0, 2 * s_r, dx);
    this->Nal = Get_Length(0, 2 * s_al, dx);
    if (this->Na % 2 == 0)
    {
        this->Na += 1;
    }
    if (this->Nr % 2 == 0)
    {
        this->Nr += 1;
    }
    if (this->Nal % 2 == 0)
    {
        this->Nal += 1;
    }
}

void Solver::set_dt(double new_dt)
{
    this->dt = new_dt;
    this->Nt = Get_Length(0, this->final_time, this->dt);
}

void Solver::set_tolerance(double new_tolerance)
{
    this->tolerance = new_tolerance;
}
