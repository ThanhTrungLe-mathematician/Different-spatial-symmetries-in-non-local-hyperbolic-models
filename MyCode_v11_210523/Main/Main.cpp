
#include "solver.h"

void run_a_test(string input, int number_of_threads);

void run_a_test_change_amplitude(string input, double amplitude, int output_position, int number_of_threads);

void run_n_test(int number_of_threads);

void run_n_test_change_amplitude(int number_of_threads);

void run_9_tests_from(double start, int number_of_threads);

void compute_convergence_rate_on_time(int number_of_threads);

void compute_convergence_rate_on_space(int number_of_threads);

void compare_CPU_time(int number_of_threads);

void create_amplitude(VECTOR &amplitude, double start, int number_of_test);

int main()
{
	string input;

	int number_of_threads = 16;
	cout << "number of threads = " << number_of_threads << endl;

	// run_a_test("input_0.txt", number_of_threads);
	// run_a_test("input_1.txt", number_of_threads);
	// run_a_test("input_2.txt", number_of_threads);
	// run_a_test("input_3.txt", number_of_threads);

	// run_n_test(number_of_threads);

	run_n_test_change_amplitude(number_of_threads);

	// compute_convergence_rate_on_time(number_of_threads);

	// compute_convergence_rate_on_space(number_of_threads);

	// compare_CPU_time(number_of_threads);

	return 0;
}

void run_n_test_change_amplitude(int number_of_threads)
{
	int start_position = 0;
	int number_of_test = 100;

	string input = "input_n_test.txt";

	Solver *solver_a_test;
	solver_a_test = new Solver(input.c_str());
	solver_a_test->read_computation_data();

	//double start = 0.001;
	double start = solver_a_test->get_initial_amplitude();
	
	delete solver_a_test;
	
	VECTOR Amplitude;
	create_amplitude(Amplitude, start, number_of_test);

	for (int i = start_position; i < number_of_test; i++)
	{
		run_a_test_change_amplitude(input, Amplitude[i], i, number_of_threads);
	}
}

void run_9_tests_from(double start, int number_of_threads)
{
	int start_position = 1;
	int number_of_test = 10;

	string input = "input_n_test.txt";
	
	VECTOR Amplitude;
	Resize(Amplitude, number_of_test);
	for (int i = 0; i < number_of_test; i++)
	{
		Amplitude[i] = start + 0.01*i;
	}

	for (int i = start_position; i < number_of_test; i++)
	{
		run_a_test_change_amplitude(input, Amplitude[i], i, number_of_threads);
	}
}

void run_a_test(string input, int number_of_threads)
{
	Solver *solver_a_test;
	solver_a_test = new Solver(input.c_str());
	cout << "1. Read input data for computation at " + get_local_time();
	solver_a_test->read_computation_data();
	solver_a_test->Solve_parallel(number_of_threads);
	delete solver_a_test;
}

void create_amplitude(VECTOR &Amplitude, double start, int number_of_test)
{
	Resize(Amplitude, number_of_test);

	if (start == 0.001)
	{
		Amplitude[0] = 0.001;
		for (int i = 1; i < number_of_test; i++)
		{
			Amplitude[i] = 0.1*i;
		}
	}
	else{
		for (int i = 0; i < number_of_test; i++)
		{
			Amplitude[i] = start + 0.1*i;
		}
	}
}

void run_a_test_change_amplitude(string input, double amplitude, int output_position, int number_of_threads)
{
	Solver *solver_a_test;
	solver_a_test = new Solver(input.c_str());
	cout << "1. Read input data for computation at " + get_local_time();
	solver_a_test->read_computation_data();

	solver_a_test->set_initial_amplitude(amplitude);

	string output = "Output_";
	string output_directory = OUTPUT_PATH + output + to_string(output_position);
	solver_a_test->set_output_directory(output_directory.c_str());

	solver_a_test->Solve_parallel(number_of_threads);
	delete solver_a_test;
}

void run_n_test(int number_of_threads)
{
	Scheme used_scheme[10] = {Upwind, MacCormack, Fractional_Step_Method, Quasi_Steady_Algorithm, Quasi_Steady_Algorithm, Quasi_Steady_Algorithm, Quasi_Steady_Algorithm, Quasi_Steady_Algorithm, Quasi_Steady_Algorithm, Quasi_Steady_Algorithm};
	SlopeLimiter used_slope[10] = {Non_limiter, Non_limiter, Non_limiter, Non_limiter, Centered_slope, Beam_Warming, Lax_Wendroof, Minmod, Superbee, MC_limiter};
	double used_dx[5];
	used_dx[0] = pow(2, -6);
	used_dx[1] = 0.01;
	used_dx[2] = pow(2, -7);
	used_dx[3] = 0.005;
	used_dx[4] = pow(2, -8);

	string input = "input_n_test.txt";
	int start_position = 0;
	int number_of_test = 22;

	double amplitude[number_of_test];
	amplitude[0] = 22.0;
	for (int i = 1; i < number_of_test; i++)
	{
		amplitude[i] = amplitude[i - 1] + 0.1;
	}

	for (int i = start_position; i < number_of_test + start_position; i++)
	{
		Solver *solver_a_test;
		solver_a_test = new Solver(input.c_str());
		cout << "1. Read input data for computation at " + get_local_time();
		solver_a_test->read_computation_data();

		solver_a_test->set_initial_amplitude(amplitude[i]);
		// solver_a_test->set_scheme(used_scheme[i]);
		// solver_a_test->set_slope_limiter(used_slope[i]);
		// solver_a_test->set_dx(used_dx[i]);
		// solver_a_test->set_dt(2.0*used_dx[i]);

		string output = "Output_";
		string output_directory = OUTPUT_PATH + output + to_string(i);
		solver_a_test->set_output_directory(output_directory.c_str());

		solver_a_test->Solve_parallel(number_of_threads);
		delete solver_a_test;
	}
}

void compute_convergence_rate_on_time(int number_of_threads)
{
	MATRIX Error;
	int number_test = 4;

	string input = "input_convergence.txt";

	Resize(Error, 2, number_test);
	string error_out = "Error/error.txt";
	error_out = OUTPUT_PATH + error_out;
	Write_Data(Error, error_out.c_str());

	MATRIX CPU_time;
	Resize(CPU_time, 2, number_test);
	string CPU_time_out = "Error/CPU_time.txt";
	CPU_time_out = OUTPUT_PATH + CPU_time_out;
	Write_Data(CPU_time, CPU_time_out.c_str());

	VECTOR result_total_density[number_test];

	for (int i = 0; i < number_test; i++)
	{
		Solver *solver_a_test;
		solver_a_test = new Solver(input.c_str());
		cout << "1. Read input data for computation at " + get_local_time();
		solver_a_test->read_computation_data();

		solver_a_test->set_dt(pow(2, i - 4 - number_test));

		string output = "Output_";
		string output_directory = OUTPUT_PATH + output + to_string(i);
		solver_a_test->set_output_directory(output_directory.c_str());

		solver_a_test->Solve_parallel(number_of_threads);

		Error[0][i] = solver_a_test->get_dt();
		result_total_density[i] = solver_a_test->get_Total_density_final_time();

		CPU_time[0][i] = solver_a_test->get_dt();
		CPU_time[1][i] = solver_a_test->get_CPU_time();

		delete solver_a_test;
	}

	for (int i = 1; i < number_test; i++)
	{
		Error[1][i] = compute_L1_discrete_error(result_total_density[i], result_total_density[0]) * 10;
	}

	Write_Data(Error, error_out.c_str());
	Write_Data(CPU_time, CPU_time_out.c_str());

	Destroy(Error);
	Destroy(CPU_time);
	for (int i = 0; i < result_total_density->size(); i++)
	{
		Destroy(result_total_density[i]);
	}
}

void compute_convergence_rate_on_space(int number_of_threads)
{
	MATRIX Error;
	int number_test = 4;

	string input = "input_convergence.txt";

	Resize(Error, 2, number_test);
	string error_out = "Error/error.txt";
	error_out = OUTPUT_PATH + error_out;
	Write_Data(Error, error_out.c_str());

	MATRIX CPU_time;
	Resize(CPU_time, 2, number_test);
	string CPU_time_out = "Error/CPU_time.txt";
	CPU_time_out = OUTPUT_PATH + CPU_time_out;
	Write_Data(CPU_time, CPU_time_out.c_str());

	VECTOR result_total_density[number_test];

	for (int i = 0; i < number_test; i++)
	{
		Solver *solver_a_test;
		solver_a_test = new Solver(input.c_str());
		cout << "1. Read input data for computation at " + get_local_time();
		solver_a_test->read_computation_data();

		solver_a_test->set_dt(pow(2, i - 4 - number_test));
		solver_a_test->set_dx(pow(2, i - 6 - number_test));

		string output = "Output_";
		string output_directory = OUTPUT_PATH + output + to_string(i);
		solver_a_test->set_output_directory(output_directory.c_str());

		solver_a_test->Solve_parallel(number_of_threads);

		Error[0][i] = solver_a_test->get_dt();
		result_total_density[i] = solver_a_test->get_Total_density_final_time();

		CPU_time[0][i] = solver_a_test->get_dt();
		CPU_time[1][i] = solver_a_test->get_CPU_time();

		delete solver_a_test;
	}

	for (int i = 1; i < number_test; i++)
	{
		Error[1][i] = compute_L1_discrete_error(result_total_density[i], result_total_density[0]) * 10;
	}

	Write_Data(Error, error_out.c_str());
	Write_Data(CPU_time, CPU_time_out.c_str());

	Destroy(Error);
	Destroy(CPU_time);
	for (int i = 0; i < result_total_density->size(); i++)
	{
		Destroy(result_total_density[i]);
	}
}

double compute_cpu_time_a_test(string input, string output, double dx, double dt, Scheme scheme, SlopeLimiter slope, int number_of_threads)
{
	Solver *solver_a_test;
	solver_a_test = new Solver(input.c_str());
	cout << "1. Read input data for computation at " + get_local_time();
	solver_a_test->read_computation_data();

	string output_directory = OUTPUT_PATH + output;
	solver_a_test->set_output_directory(output_directory.c_str());

	solver_a_test->set_dx(dx);
	solver_a_test->set_dt(dt);
	solver_a_test->set_scheme(scheme);
	solver_a_test->set_slope_limiter(slope);

	solver_a_test->Solve_parallel(number_of_threads);
	double cpu_time = solver_a_test->get_CPU_time();
	delete solver_a_test;

	return cpu_time;
}

void compare_CPU_time(int number_of_threads)
{
	MATRIX Error;

	string input = "input_CPU_time.txt";

	MATRIX CPU_time;
	Resize(CPU_time, 7, 3);
	string CPU_time_out = "Error/CPU_time.txt";
	CPU_time_out = OUTPUT_PATH + CPU_time_out;
	Write_Data(CPU_time, CPU_time_out.c_str());

	double dx[3] = {pow(2, -7), pow(2, -7), pow(2, -7)};
	double dt[3] = {pow(2, -5), pow(2, -6), pow(2, -7)};

	CPU_time[0][0] = dt[0];
	CPU_time[0][1] = dt[1];
	CPU_time[0][2] = dt[2];

	CPU_time[1][0] = compute_cpu_time_a_test(input, "Output_0", dx[0], dt[0], Upwind, Non_limiter, number_of_threads);
	CPU_time[1][1] = compute_cpu_time_a_test(input, "Output_1", dx[1], dt[1], Upwind, Non_limiter, number_of_threads);
	CPU_time[1][2] = compute_cpu_time_a_test(input, "Output_2", dx[2], dt[2], Upwind, Non_limiter, number_of_threads);

	CPU_time[2][0] = compute_cpu_time_a_test(input, "Output_3", dx[0], dt[0], MacCormack, Non_limiter, number_of_threads);
	CPU_time[2][1] = compute_cpu_time_a_test(input, "Output_4", dx[1], dt[1], MacCormack, Non_limiter, number_of_threads);
	CPU_time[2][2] = compute_cpu_time_a_test(input, "Output_5", dx[2], dt[2], MacCormack, Non_limiter, number_of_threads);

	CPU_time[3][0] = compute_cpu_time_a_test(input, "Output_6", dx[0], dt[0], Fractional_Step_Method, Non_limiter, number_of_threads);
	CPU_time[3][1] = compute_cpu_time_a_test(input, "Output_7", dx[1], dt[1], Fractional_Step_Method, Non_limiter, number_of_threads);
	CPU_time[3][2] = compute_cpu_time_a_test(input, "Output_8", dx[2], dt[2], Fractional_Step_Method, Non_limiter, number_of_threads);

	CPU_time[4][0] = compute_cpu_time_a_test(input, "Output_9", dx[0], dt[0], Quasi_Steady_Algorithm, Non_limiter, number_of_threads);
	CPU_time[4][1] = compute_cpu_time_a_test(input, "Output_10", dx[1], dt[1], Quasi_Steady_Algorithm, Non_limiter, number_of_threads);
	CPU_time[4][2] = compute_cpu_time_a_test(input, "Output_11", dx[2], dt[2], Quasi_Steady_Algorithm, Non_limiter, number_of_threads);

	CPU_time[5][0] = compute_cpu_time_a_test(input, "Output_12", dx[0], dt[0], Quasi_Steady_Algorithm, Beam_Warming, number_of_threads);
	CPU_time[5][1] = compute_cpu_time_a_test(input, "Output_13", dx[1], dt[1], Quasi_Steady_Algorithm, Beam_Warming, number_of_threads);
	CPU_time[5][2] = compute_cpu_time_a_test(input, "Output_14", dx[2], dt[2], Quasi_Steady_Algorithm, Beam_Warming, number_of_threads);

	CPU_time[6][0] = compute_cpu_time_a_test(input, "Output_15", dx[0], dt[0], Quasi_Steady_Algorithm, MC_limiter, number_of_threads);
	CPU_time[6][1] = compute_cpu_time_a_test(input, "Output_16", dx[1], dt[1], Quasi_Steady_Algorithm, MC_limiter, number_of_threads);
	CPU_time[6][2] = compute_cpu_time_a_test(input, "Output_17", dx[2], dt[2], Quasi_Steady_Algorithm, MC_limiter, number_of_threads);

	Write_Data(CPU_time, CPU_time_out.c_str());

	Destroy(CPU_time);
}