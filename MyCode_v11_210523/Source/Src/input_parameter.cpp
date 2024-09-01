//
//This file contains all methods of class Input_parameter
//

#include "input_parameter.h"

Input_parameter::Input_parameter(){}

Input_parameter::~Input_parameter(){
    data.clear();
}

Input_parameter::Input_parameter(string filename){
    filename = INPUT_PATH + filename;
    this->Parser_textfile(filename.c_str());

    this->left_boundary = atof(this->Get_Value_of_variable("<left_boundary>").c_str());

    this->right_boundary = atof(this->Get_Value_of_variable("<right_boundary>").c_str());

    this->dx = atof(this->Get_Value_of_variable("<dx>").c_str());

    this->final_time = atof(this->Get_Value_of_variable("<final_time>").c_str());

    this->dt = atof(this->Get_Value_of_variable("<dt>").c_str());

    this->output_directory = (OUTPUT_PATH + this->Get_Value_of_variable("<output_directory>")).c_str();

    this->gamma = atof(this->Get_Value_of_variable("<gamma>").c_str());

    this->lambda_1 = atof(this->Get_Value_of_variable("<lambda_1>").c_str());

    this->lambda_2 = atof(this->Get_Value_of_variable("<lambda_2>").c_str());

    this->q_r = atof(this->Get_Value_of_variable("<q_r>").c_str());

    this->q_a = atof(this->Get_Value_of_variable("<q_a>").c_str());

    this->q_al = atof(this->Get_Value_of_variable("<q_al>").c_str());

    this->initial_amplitude = atof(this->Get_Value_of_variable("<initial_amplitude>").c_str());

    this->model = (this->Get_Value_of_variable("<model>")).c_str();

    this->scheme = (this->Get_Value_of_variable("<scheme>")).c_str();

    this->slope_limiter = (this->Get_Value_of_variable("<slope_limiter>")).c_str();

    this->boundary_condition = (this->Get_Value_of_variable("<boundary_condition>")).c_str();

    this->tolerance = atof(this->Get_Value_of_variable("<tolerance>").c_str());
}

string Input_parameter::get_output_directory(){
    return this->output_directory;
}

double Input_parameter::get_left_boundary(){
    return this->left_boundary;
}

double Input_parameter::get_right_boundary(){
    return this->right_boundary;
}

double Input_parameter::get_dx(){
    return this->dx;
}

double Input_parameter::get_final_time(){
    return this->final_time;
}

double Input_parameter::get_dt(){
    return this->dt;
}

double Input_parameter::get_gamma(){
    return this->gamma;
}

double Input_parameter::get_lambda_1(){
    return this->lambda_1;
}

double Input_parameter::get_lambda_2(){
    return this->lambda_2;
}

double Input_parameter::get_q_r(){
    return this->q_r;
}

double Input_parameter::get_q_a(){
    return this->q_a;
}

double Input_parameter::get_q_al(){
    return this->q_al;
}

double Input_parameter::get_initial_amplitude(){
    return this->initial_amplitude;
}

string Input_parameter::get_model(){
    return this->model;
}

string Input_parameter::get_scheme(){
    return this->scheme;
}

string Input_parameter::get_slope_limiter(){
    return this->slope_limiter;
}

string Input_parameter::get_boundary_condition(){
    return this->boundary_condition;
}

double Input_parameter::get_tolerance(){
    return this->tolerance;
}


void Input_parameter::Parser_textfile(string filename){
    ifstream readfile(filename, ios::in);

    //if the file can not be open then stop
    if(!readfile){
        cerr<<"Impossible to open the "<<filename<<" file!"<<endl;
        exit(EXIT_FAILURE);
    }

    this->number_of_lines = 0;
    string teststring;
    //counting the number of lines in the file
	while(!readfile.eof()){
        getline(readfile, teststring);
		this->number_of_lines ++;
    }

    readfile.clear();
    readfile.seekg(0, ios::beg); //the file will be read again from the beginning

    this->data.resize(this->number_of_lines);

    // data is a tabular where the useful lines of the input parameter are copied
	// the comments (with a #) and the empty lines are not copied 
    int j = -1;
    int i = 0;

	while (!readfile.eof()){
		getline (readfile, teststring);		// the whole line is copied in the teststring
		j = int(teststring.find("#"));      // j is the position of the first #

		if (teststring!="" && j!=0){	    // if the line is not empty and if it is not a comment line (which begins with a #)
			if (j>0){					    // if there is a comment (which begins with a #) after the datas
				teststring.erase (teststring.begin() + j, teststring.end());	// erase the comment
			}
			this->data[i] = teststring;
			i++;
		}
	}
	
	this->number_of_lines = i;					// number of non-empty lines of data
    
	readfile.close();
    //cout<<"The input parameters in file "<<filename<<" has been read!"<<endl;
}

string Input_parameter::Get_Value_of_variable(string name_of_variable){
    int position_line = -1; //the position of the line that contains the name of variable
    int first_position = -1;

    // file the position of the line that contains the name of variable
    for(int i = 0; i < this->number_of_lines; i++){
        first_position = int(data[i].find(name_of_variable)); //int(find(...)) return -1 if can not find the value
        if(first_position > 0){
            position_line = i;
            break;
        }
    }

    //If the string name_of_variable is not in the parameters files then stop
    if(first_position <= -1 && position_line <= -1){
        cerr<<"No entry for the variable "<<name_of_variable<<endl;
        exit(EXIT_FAILURE);
    }

    string temp = data[position_line];

    first_position = int(temp.find("::"));
    if (first_position <= -1){					// if the :: are not found
		cerr << "Bad syntax for "<<name_of_variable<<". The syntax is: description <variable>:: value"<<endl;
		exit(EXIT_FAILURE);
	}

    
    // erase the description of the variable and ::
    temp.erase(0, first_position+2);	
	
    // erase the white spaces before the string
	temp.erase(0, int(temp.find_first_not_of(" ")));		
    
    // erase the white spaces after the string
	temp.erase(temp.begin() + int(temp.find_last_not_of(" ")) + 1, temp.end());	

    return temp;
}