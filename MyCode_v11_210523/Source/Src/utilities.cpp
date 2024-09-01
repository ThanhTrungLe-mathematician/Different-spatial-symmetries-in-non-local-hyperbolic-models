//
//	This file contains all global functions
//

#include "utilities.h"

void Randomize()
{
    srand((unsigned)time(NULL));
}

double get_cpu_time_sec(){
    return (double) clock () / (double) CLOCKS_PER_SEC;
}

double get_cpu_time_hour(){
    return (double) clock () / (double) CLOCKS_PER_SEC / 3600.0;
}

string get_local_time(){
    time_t now = time(0);
    string local_time = ctime(&now);
    return local_time;
}

double Rand()
{
    return (double)rand() / (RAND_MAX);
}

int Get_Length(double start_point, double end_point, double step)
{
    if (end_point <= start_point || step <= 0)
    {
        cerr << "The parameters: start_point, end_point, step are not correct!" << endl;
        exit(EXIT_FAILURE);
    }
    /*double temp = start_point;
    int length = 1;
    while(temp < end_point){
        cout<<temp<<endl;
        temp += step;
        length += 1;
        cout<<temp<<"\t"<<length<<endl;
    }
    return length;*/
    return int((end_point - start_point) / step + 1);
}

void Create_Vector(VECTOR &aVector, double start_point, double end_point, double step)
{
    int length = Get_Length(start_point, end_point, step);
    Resize(aVector, length);
    for (int i = 0; i < length; ++i)
    {
        aVector[i] = start_point + i * step;
    }
}

void Resize(VECTOR &aVector, int Nx)
{
    if (Nx <= 0)
    {
        cerr << "The size of Vector is " << Nx << ". This must be positive!" << endl;
        exit(EXIT_FAILURE);
    }

    aVector.resize(Nx);
}

void Resize(MATRIX &aMatrix, int Nx, int Ny)
{
    if (Nx <= 0)
    {
        cerr << "The first parameter in the size of Matrix is " << Nx << ". This must be positive!" << endl;
        exit(EXIT_FAILURE);
    }
    if (Ny <= 0)
    {
        cerr << "The second parameter in the size of Matrix is " << Ny << ". This must be positive!" << endl;
        exit(EXIT_FAILURE);
    }

    aMatrix.resize(Nx);
    for (int i = 0; i < Nx; ++i)
    {
        Resize(aMatrix[i], Ny);
    }
}

void Resize(TENSOR3D &aTensor3D, int Nt, int Nx, int Ny)
{
    if (Nt <= 0)
    {
        cerr << "The first parameter in the size of Tensor3D is " << Nt << ". This must be positive!" << endl;
        exit(EXIT_FAILURE);
    }
    if (Nx <= 0)
    {
        cerr << "The second parameter in the size of Tensor3D is " << Nx << ". This must be positive!" << endl;
        exit(EXIT_FAILURE);
    }
    if (Ny <= 0)
    {
        cerr << "The third parameter in the size of Tensor3D is " << Ny << ". This must be positive!" << endl;
        exit(EXIT_FAILURE);
    }

    aTensor3D.resize(Nt);
    for (int i = 0; i < Nt; ++i)
    {
        Resize(aTensor3D[i], Nx, Ny);
    }
}

void Resize(TENSOR4D &aTensor4D, int Nt, int Na, int Nx, int Ny)
{
    if (Nt <= 0)
    {
        cerr << "The first parameter in the size of Tensor4D is " << Nt << ". This must be positive!" << endl;
        exit(EXIT_FAILURE);
    }
    if (Nx <= 0)
    {
        cerr << "The second parameter in the size of Tensor4D is " << Na << ". This must be positive!" << endl;
        exit(EXIT_FAILURE);
    }
    if (Ny <= 0)
    {
        cerr << "The third parameter in the size of Tensor4D is " << Nx << ". This must be positive!" << endl;
        exit(EXIT_FAILURE);
    }

    if (Ny <= 0)
    {
        cerr << "The fourth parameter in the size of Tensor4D is " << Ny << ". This must be positive!" << endl;
        exit(EXIT_FAILURE);
    }

    aTensor4D.resize(Nt);
    for (int i = 0; i < Nt; ++i)
    {
        Resize(aTensor4D[i], Na, Nx, Ny);
    }
}

void Destroy(VECTOR &aVector)
{
    aVector.clear();
}

void Destroy(MATRIX &aMatrix)
{
    for (size_t i = 0; i < aMatrix.size(); ++i)
    {
        Destroy(aMatrix[i]);
    }
    aMatrix.clear();
}

void Destroy(TENSOR3D &aTensor3D)
{
    for (size_t i = 0; i < aTensor3D.size(); ++i)
    {
        Destroy(aTensor3D[i]);
    }
    aTensor3D.clear();
}

void Destroy(TENSOR4D &aTensor4D)
{
    for (size_t i = 0; i < aTensor4D.size(); ++i)
    {
        Destroy(aTensor4D[i]);
    }
    aTensor4D.clear();
}

void Print(VECTOR aVector)
{
    int length = aVector.size();
    for (int i = 0; i < length; ++i)
    {
        cout << aVector[i];
        if (i < length - 1)
        {
            cout << "\t";
        }
    }
}

void Print(MATRIX aMatrix)
{
    int row = aMatrix.size();
    for (int i = 0; i < row; ++i)
    {
        Print(aMatrix[i]);
        if (i < row - 1)
        {
            cout << endl;
        }
    }
}

void Print(TENSOR3D aTensor3D)
{
    int height = aTensor3D.size();
    for (int i = 0; i < height; ++i)
    {
        cout << "#" << endl;
        Print(aTensor3D[i]);
        cout << endl;
    }
    cout << "#";
}

void Write_Data(VECTOR &aVector, string filename)
{
    ofstream writefile(filename, ios::out);

    // if the file can not be open then stop
    if (!writefile)
    {
        cerr << "Impossible to open the " << filename << " file!" << endl;
        exit(EXIT_FAILURE);
    }

    writefile.precision(PRECISION); // set precision
    int length = aVector.size();
    for (int i = 0; i < length; ++i)
    {
        writefile << aVector[i];
        if (i < length - 1)
        {
            writefile << "\t";
        }
    }
    writefile.close();
    //cout << "The Vector has been written to: " << filename << endl;
}

void Write_Data(VECTOR &aVector, VECTOR &bVector, string filename)
{
    if(aVector.size() != bVector.size())
    {
        cerr<<"The size of two vectors must be equal! Please check again!"<<endl;
        exit(EXIT_FAILURE);
    }

    ofstream writefile(filename, ios::out);

    // if the file can not be open then stop
    if (!writefile)
    {
        cerr << "Impossible to open the " << filename << " file!" << endl;
        exit(EXIT_FAILURE);
    }

    writefile.precision(PRECISION); // set precision
    int length = aVector.size();
    for (int i = 0; i < length; ++i)
    {
        writefile << aVector[i] + bVector[i];
        if (i < length - 1)
        {
            writefile << "\t";
        }
    }
    writefile.close();
    //cout << "The Vector has been written to: " << filename << endl;
}

void Write_Data(MATRIX &aMatrix, string filename)
{
    ofstream writefile(filename, ios::out);

    // if the file can not be open then stop
    if (!writefile)
    {
        cerr << "Impossible to open the " << filename << " file!" << endl;
        exit(EXIT_FAILURE);
    }

    writefile.precision(PRECISION); // set precision
    int row = aMatrix.size();
    int column = aMatrix[0].size();

    for (int i = 0; i < row; ++i)
    {
        for (int j = 0; j < column; ++j)
        {
            writefile << aMatrix[i][j];
            if (j < column - 1)
            {
                writefile << "\t";
            }
        }
        if (i < row - 1)
        {
            writefile << endl;
        }
    }
    writefile.close();
    //cout << "The Matrix has been written to: " << filename << endl;
}

void Write_Data(MATRIX &aMatrix, MATRIX &bMatrix, string filename)
{
    if((aMatrix.size() != bMatrix.size()) || (aMatrix[0].size() != bMatrix[0].size()))
    {
        cerr<<"The size of two matrixs must be equal! Please check again!"<<endl;
        exit(EXIT_FAILURE);
    }

    ofstream writefile(filename, ios::out);

    // if the file can not be open then stop
    if (!writefile)
    {
        cerr << "Impossible to open the " << filename << " file!" << endl;
        exit(EXIT_FAILURE);
    }

    writefile.precision(PRECISION); // set precision
    int row = aMatrix.size();
    int column = aMatrix[0].size();

    for (int i = 0; i < row; ++i)
    {
        for (int j = 0; j < column; ++j)
        {
            writefile << aMatrix[i][j] + bMatrix[i][j];
            if (j < column - 1)
            {
                writefile << "\t";
            }
        }
        if (i < row - 1)
        {
            writefile << endl;
        }
    }
    writefile.close();
    //cout << "The Matrix has been written to: " << filename << endl;
}

void Write_Data(MATRIX &aMatrix, string filename, int number_of_final_step)
{
    ofstream writefile(filename, ios::out);

    // if the file can not be open then stop
    if (!writefile)
    {
        cerr << "Impossible to open the " << filename << " file!" << endl;
        exit(EXIT_FAILURE);
    }

    writefile.precision(PRECISION); // set precision
    int row = aMatrix.size();
    int column = aMatrix[0].size();

    for (int i = row - number_of_final_step - 1; i < row; ++i)
    {
        for (int j = 0; j < column; ++j)
        {
            writefile << aMatrix[i][j];
            if (j < column - 1)
            {
                writefile << "\t";
            }
        }
        if (i < row - 1)
        {
            writefile << endl;
        }
    }
    writefile.close();
    //cout << "The Matrix has been written to: " << filename << endl;
}

void Write_Data(MATRIX &aMatrix, MATRIX &bMatrix, string filename, int number_of_final_step)
{
    if((aMatrix.size() != bMatrix.size()) || (aMatrix[0].size() != bMatrix[0].size()))
    {
        cerr<<"The size of two matrixs must be equal! Please check again!"<<endl;
        exit(EXIT_FAILURE);
    }

    ofstream writefile(filename, ios::out);

    // if the file can not be open then stop
    if (!writefile)
    {
        cerr << "Impossible to open the " << filename << " file!" << endl;
        exit(EXIT_FAILURE);
    }

    writefile.precision(PRECISION); // set precision
    int row = aMatrix.size();
    int column = aMatrix[0].size();

    for (int i = row - number_of_final_step - 1; i < row; ++i)
    {
        for (int j = 0; j < column; ++j)
        {
            writefile << aMatrix[i][j] + bMatrix[i][j];
            if (j < column - 1)
            {
                writefile << "\t";
            }
        }
        if (i < row - 1)
        {
            writefile << endl;
        }
    }
    writefile.close();
    //cout << "The Matrix has been written to: " << filename << endl;
}

void Write_Data(MATRIX &aMatrix, string filename, int scale_x, int scale_y)
{
    ofstream writefile(filename, ios::out);

    // if the file can not be open then stop
    if (!writefile)
    {
        cerr << "Impossible to open the " << filename << " file!" << endl;
        exit(EXIT_FAILURE);
    }

    writefile.precision(PRECISION); // set precision
    int row = aMatrix.size();
    int column = aMatrix[0].size();

    for (int i = 0; i < row; ++i)
    {
        if (i % scale_x == 0)
        {
            for (int j = 0; j < column; ++j)
            {
                if (j % scale_y == 0)
                {
                    writefile << aMatrix[i][j];
                    if (j < column - 1)
                    {
                        writefile << "\t";
                    }
                }
            }
            if (i < row - 1)
            {
                writefile << endl;
            }
        }
    }
    writefile.close();
    //cout << "The Matrix has been written to: " << filename << endl;
}

void Write_Data(MATRIX &aMatrix, MATRIX &bMatrix, string filename, int scale_x, int scale_y)
{
    if((aMatrix.size() != bMatrix.size()) || (aMatrix[0].size() != bMatrix[0].size()))
    {
        cerr<<"The size of two matrixs must be equal! Please check again!"<<endl;
        exit(EXIT_FAILURE);
    }

    ofstream writefile(filename, ios::out);

    // if the file can not be open then stop
    if (!writefile)
    {
        cerr << "Impossible to open the " << filename << " file!" << endl;
        exit(EXIT_FAILURE);
    }

    writefile.precision(PRECISION); // set precision
    int row = aMatrix.size();
    int column = aMatrix[0].size();

    for (int i = 0; i < row; ++i)
    {
        if (i % scale_x == 0)
        {
            for (int j = 0; j < column; ++j)
            {
                if (j % scale_y == 0)
                {
                    writefile << aMatrix[i][j] + bMatrix[i][j];
                    if (j < column - 1)
                    {
                        writefile << "\t";
                    }
                }
            }
            if (i < row - 1)
            {
                writefile << endl;
            }
        }
    }
    writefile.close();
    //cout << "The Matrix has been written to: " << filename << endl;
}

void Write_Data(TENSOR3D &aTensor3D, string filename)
{
    ofstream writefile(filename, ios::out);

    // if the file can not be open then stop
    if (!writefile)
    {
        cerr << "Impossible to open the " << filename << " file!" << endl;
        exit(EXIT_FAILURE);
    }

    writefile.precision(PRECISION); // set precision
    int height = aTensor3D.size();
    int row = aTensor3D[0].size();
    int column = aTensor3D[0][0].size();

    for (int k = 0; k < height; ++k)
    {
        writefile << "#" << endl;
        for (int i = 0; i < row; ++i)
        {
            for (int j = 0; j < column; ++j)
            {
                writefile << aTensor3D[k][i][j];
                if (j < column - 1)
                {
                    writefile << "\t";
                }
            }
            if (i < row - 1)
            {
                writefile << endl;
            }
        }
        writefile << endl;
    }
    writefile << "#";
    writefile.close();
    //cout << "The Tensor3D has been written to: " << filename << endl;
}

void Read_Data(VECTOR &aVector, string filename)
{
    ifstream readfile(filename, ios::in);

    // if the file can not be open then stop
    if (!readfile)
    {
        cerr << "Impossible to open the " << filename << " file!" << endl;
        exit(EXIT_FAILURE);
    }

    int number_of_element = 0;
    double test = 0;

    // counting the element of the file
    while (!readfile.eof())
    {
        readfile >> test;
        number_of_element += 1;
    }

    int length = aVector.size();

    if (length == 0)
    {
        Resize(aVector, number_of_element); // change the size of the Vector
    }
    else if (number_of_element != length)
    {
        cerr << "The size of vector in the file " << filename << " is: " << number_of_element << endl;
        cerr << "The size of object is: " << length << endl;
        cerr << "The two above values must be equal!" << endl;
        readfile.close();
        exit(EXIT_FAILURE);
    }

    readfile.clear();
    readfile.seekg(0, ios::beg); // the file will be read again from the beginning
    for (size_t i = 0; i < aVector.size(); ++i)
    {
        readfile >> aVector[i];
    }
    readfile.close();
    //cout << "The vector in file " << filename << " has been read!" << endl;
}

void Read_Data(MATRIX &aMatrix, string filename)
{
    ifstream readfile(filename, ios::in);

    // if the file can not be open then stop
    if (!readfile)
    {
        cerr << "Impossible to open the " << filename << " file!" << endl;
        exit(EXIT_FAILURE);
    }

    int number_of_element = 0;
    int number_of_row_file = 0;
    int number_of_column_file = 0;
    double test = 0;
    string teststring;

    // counting the element of the file
    while (!readfile.eof())
    {
        readfile >> test;
        number_of_element++;
    }

    readfile.clear();
    readfile.seekg(0, ios::beg); // the file will be read again from the beginning

    // counting number of row in the file
    while (!readfile.eof())
    {
        getline(readfile, teststring);
        number_of_row_file += 1;
    }

    number_of_column_file = number_of_element / number_of_row_file;

    int NumberOfRow = aMatrix.size();
    int NumberOfColumn = 0;
    if (NumberOfRow != 0)
    {
        NumberOfColumn = aMatrix[0].size();
    }

    if (NumberOfRow == 0 && NumberOfColumn == 0)
    {
        Resize(aMatrix, number_of_row_file, number_of_column_file); // change the size of the Matrix
    }
    else if (number_of_row_file != NumberOfRow)
    {
        cerr << "The first parameter in the size of Matrix in the file " << filename << " is: " << number_of_row_file << endl;
        cerr << "The first parameter in the size of Matrix object is: " << NumberOfRow << endl;
        cerr << "The two above values must be equal!" << endl;
        readfile.close();
        exit(EXIT_FAILURE);
    }
    else if (number_of_column_file != NumberOfColumn)
    {
        cerr << "The second parameter in the size of Matrix in the file " << filename << " is: " << number_of_column_file << endl;
        cerr << "The second parameter in the size of Matrix object is: " << NumberOfColumn << endl;
        cerr << "The two above values must be equal!" << endl;
        readfile.close();
        exit(EXIT_FAILURE);
    }

    readfile.clear();
    readfile.seekg(0, ios::beg); // the file will be read again from the beginning

    for (size_t i = 0; i < aMatrix.size(); ++i)
    {
        for (size_t j = 0; j < aMatrix[0].size(); ++j)
        {
            readfile >> aMatrix[i][j];
        }
    }
    readfile.close();
    //cout << "The matrix in file " << filename << " has been read!" << endl;
}

void Read_Data(TENSOR3D &aTensor3D, string filename)
{
    ifstream readfile(filename, ios::in);

    // if the file can not be open then stop
    if (!readfile)
    {
        cerr << "Impossible to open the " << filename << " file!" << endl;
        exit(EXIT_FAILURE);
    }

    int number_of_element = 0;
    string teststring;

    // counting the number of element of one block matrix of Tensor3D in the file
    while (!readfile.eof())
    {
        readfile >> teststring;
        number_of_element += 1;
        if (number_of_element == 1 && teststring != "#")
        {
            cerr << "The structure of the Tensor3D in the " << filename << " file is not correct!" << endl;
            cerr << "It must be start with #!" << endl;
            readfile.close();
            exit(EXIT_FAILURE);
        }
        if (number_of_element > 1 && teststring == "#")
        {
            break;
        }
    }

    number_of_element -= 2;

    readfile.clear();
    readfile.seekg(0, ios::beg); // the file will be read again from the beginning

    // counting the length of block Matrix and length of Tensor3D in the file
    int length_of_block = 0;
    int number_of_block = 0;
    while (!readfile.eof())
    {
        getline(readfile, teststring);
        if (teststring == "#")
        {
            number_of_block += 1;
        }
        if (number_of_block < 2)
        {
            length_of_block += 1;
        }
    }

    length_of_block -= 1;
    number_of_block -= 1;

    // compute the size (first_para, second_para, third_para) of Tensor3D in the file

    int first_para_file = number_of_block;
    int second_para_file = length_of_block;
    int third_para_file = number_of_element / length_of_block;

    // compute the size of Tensor3D object
    int first_size = aTensor3D.size();
    int second_size = 0;
    int third_size = 0;
    if (first_size != 0)
    {
        second_size = aTensor3D[0].size();
        third_size = aTensor3D[0][0].size();
    }

    if (first_size == 0 && second_size == 0 && third_size == 0)
    {
        Resize(aTensor3D, first_para_file, second_para_file, third_para_file);
    }
    else if (first_size != first_para_file)
    {
        cerr << "The first parameter in the size of Tensor3D in the file " << filename << " is: " << first_para_file << endl;
        cerr << "The first parameter in the size of Tensor3D object is: " << first_size << endl;
        cerr << "The two above values must be equal!" << endl;
        readfile.close();
        exit(EXIT_FAILURE);
    }
    else if (second_size != second_para_file)
    {
        cerr << "The second parameter in the size of Tensor3D in the file " << filename << " is: " << second_para_file << endl;
        cerr << "The second parameter in the size of Tensor3D object is: " << second_size << endl;
        cerr << "The two above values must be equal!" << endl;
        readfile.close();
        exit(EXIT_FAILURE);
    }
    else if (third_size != third_para_file)
    {
        cerr << "The third parameter in the size of Tensor3D in the file " << filename << " is: " << third_para_file << endl;
        cerr << "The third parameter in the size of Tensor3D object is: " << third_size << endl;
        cerr << "The two above values must be equal!" << endl;
        readfile.close();
        exit(EXIT_FAILURE);
    }

    readfile.clear();
    readfile.seekg(0, ios::beg); // the file will be read again from the beginning

    for (size_t i = 0; i < aTensor3D.size(); ++i)
    {
        readfile >> teststring; // read "#"
        for (size_t j = 0; j < aTensor3D[0].size(); ++j)
        {
            for (size_t k = 0; k < aTensor3D[0][0].size(); ++k)
            {
                readfile >> aTensor3D[i][j][k];
            }
        }
    }
    readfile >> teststring; // read "#"
    readfile.close();
    //cout << "The Tensor3D in file " << filename << " has been read!" << endl;
}

double Get_Max_Value(VECTOR &aVector)
{
    if (aVector.size() <= 0)
    {
        cerr << "The size of Vector is zero, so the maximum value can not be found!" << endl;
        exit(EXIT_FAILURE);
    }
    return *max_element(aVector.begin(), aVector.end());
}

double Get_Max_Value(MATRIX &aMatrix)
{
    if (aMatrix.size() <= 0)
    {
        cerr << "The size of Matrix is zero, so the maximum value can not be found!" << endl;
        exit(EXIT_FAILURE);
    }

    VECTOR temp;
    Resize(temp, aMatrix.size());
    for (size_t i = 0; i < temp.size(); ++i)
    {
        temp[i] = Get_Max_Value(aMatrix[i]);
    }
    double max_value = Get_Max_Value(temp);
    Destroy(temp);
    return max_value;
}

double Get_Max_Value(TENSOR3D &aTensor3D)
{
    if (aTensor3D.size() <= 0)
    {
        cerr << "The size of Tensor3D is zero, so the maximum value can not be found!" << endl;
        exit(EXIT_FAILURE);
    }
    VECTOR temp;
    Resize(temp, aTensor3D.size());
    for (size_t i = 0; i < temp.size(); ++i)
    {
        temp[i] = Get_Max_Value(aTensor3D[i]);
    }
    double max_value = Get_Max_Value(temp);
    Destroy(temp);
    return max_value;
}

double Get_Max_Value(TENSOR4D &aTensor4D)
{
    if (aTensor4D.size() <= 0)
    {
        cerr << "The size of Tensor4D is zero, so the maximum value can not be found!" << endl;
        exit(EXIT_FAILURE);
    }
    VECTOR temp;
    Resize(temp, aTensor4D.size());
    for (size_t i = 0; i < temp.size(); ++i)
    {
        temp[i] = Get_Max_Value(aTensor4D[i]);
    }
    double max_value = Get_Max_Value(temp);
    Destroy(temp);
    return max_value;
}

void Create_Vector(VECTOR &aVector, double (*function)(double), VECTOR V_x)
{
    int length = V_x.size();
    Resize(aVector, length);
    for (int i = 0; i < length; ++i)
    {
        aVector[i] = function(V_x[i]);
    }
}

void Create_Matrix(MATRIX &aMatrix, double (*function)(double, double), VECTOR V_x, VECTOR V_y)
{
    int length_Vx = V_x.size();
    int length_Vy = V_y.size();
    Resize(aMatrix, length_Vx, length_Vy);
    for (int i = 0; i < length_Vx; ++i)
    {
        for (int j = 0; j < length_Vy; ++j)
        {
            aMatrix[i][j] = function(V_x[i], V_y[j]);
        }
    }
}

void Create_Tensor3D(TENSOR3D &aTensor3D, double (*function)(double, double, double), VECTOR V_t, VECTOR V_x, VECTOR V_y)
{
    int length_Vt = V_t.size();
    int length_Vx = V_x.size();
    int length_Vy = V_y.size();
    Resize(aTensor3D, length_Vt, length_Vx, length_Vy);

    for (int i = 0; i < length_Vt; ++i)
    {
        for (int j = 0; j < length_Vx; ++j)
        {
            for (int k = 0; k < length_Vy; ++k)
            {
                aTensor3D[i][j][k] = function(V_t[i], V_x[j], V_y[k]);
            }
        }
    }
}

void Create_Tensor4D(TENSOR4D &aTensor4D, double (*function)(double, double, double, double), VECTOR V_t, VECTOR V_a, VECTOR V_x, VECTOR V_y)
{
    int length_Vt = V_t.size();
    int length_Va = V_a.size();
    int length_Vx = V_x.size();
    int length_Vy = V_y.size();
    Resize(aTensor4D, length_Vt, length_Va, length_Vx, length_Vy);

    for (int i = 0; i < length_Vt; ++i)
    {
        for (int j = 0; j < length_Va; ++j)
        {
            for (int k = 0; k < length_Vx; ++k)
            {
                for (int l = 0; l < length_Vy; ++l)
                {
                    aTensor4D[i][j][k][l] = function(V_t[i], V_a[j], V_x[k], V_y[l]);
                }
            }
        }
    }
}

double compute_L1_discrete_error(TENSOR3D &Numerical_sol, TENSOR3D &Reference_sol)
{
    double sum_error = 0;
    double sum_norm = 0;

    int scale_h = (Reference_sol.size() - 1) / (Numerical_sol.size() - 1);
    int scale_i = (Reference_sol[0].size() - 1) / (Numerical_sol[0].size() - 1);
    int scale_j = (Reference_sol[0][0].size() - 1) / (Numerical_sol[0][0].size() - 1);

    for (int h = 0; h < Numerical_sol.size(); ++h)
    {
        for (int i = 0; i < Numerical_sol[0].size(); ++i)
        {
            for (int j = 0; j < Numerical_sol[0][0].size(); ++j)
            {
                sum_error += abs(Numerical_sol[h][i][j] - Reference_sol[h * scale_h][i * scale_i][j * scale_j]);
                sum_norm += Reference_sol[h * scale_h][i * scale_i][j * scale_j];
            }
        }
    }

    return sum_error / sum_norm;
}

double compute_L1_discrete_error(MATRIX &Numerical_sol, MATRIX &Reference_sol)
{

    double sum_error = 0;
    double sum_norm = 0;

    int scale_i = (Reference_sol.size() - 1) / (Numerical_sol.size() - 1);
    int scale_j = (Reference_sol[0].size() - 1) / (Numerical_sol[0].size() - 1);

    // cout<<scale_i<<"\t"<<scale_j<<endl;

    for (int i = 0; i < Numerical_sol.size(); ++i)
    {
        for (int j = 0; j < Numerical_sol[0].size(); ++j)
        {
            sum_error += abs(Numerical_sol[i][j] - Reference_sol[i * scale_i][j * scale_j]);
            sum_norm += Reference_sol[i * scale_i][j * scale_j];

            // if(abs(Numerical_sol[i][j] - Reference_sol[i*scale_i][j*scale_j])>1){
            //         cout<<i<<"\t"<<j<<"\t"<<Numerical_sol[i][j]<<"\t"<<Reference_sol[i*scale_i][j*scale_j]<<endl;
            //     }
        }
        // cout<<i<<"\t"<<i*scale_i<<endl;
    }

    // cout<<sum_error<<"\t"<<sum_norm<<endl;
    return sum_error / sum_norm;
}

double compute_L1_discrete_error(VECTOR &First_vector, VECTOR &Second_vector)
{
    double sum_error = 0;
    int length = First_vector.size();

    if(length != Second_vector.size())
    {
        cerr<<"The length of the first vector is: "<<First_vector.size()<<endl;
        cerr<<"The length of the second vector is: "<<Second_vector.size()<<endl;
        cerr<<"The length of two vectors must be equal!"<<endl;
        exit(EXIT_FAILURE);
    }

    for(int i = 0; i <length; i++)
    {
        sum_error += abs(First_vector[i] - Second_vector[i]);
    }
    return sum_error/length;
}

double compute_L1_discrete_error(VECTOR &A_vector, VECTOR &B_vector, VECTOR &C_vector, VECTOR &D_vector)
{
    double sum_error = 0;
    int length = A_vector.size();

    for(int i = 0; i <length; i++)
    {
        sum_error += abs(A_vector[i] + B_vector[i] - C_vector[i] - D_vector[i]);
    }
    return sum_error/length;
}

double min_of(double a, double b)
{
    if(a < b){
        return a;
    }
    else{
        return b;
    }
}

double min_of(double a, double b, double c)
{
    return min_of(min_of(a,b), c);
}

double max_of(double a, double b)
{
    if(a > b)
    {
        return a;
    }
    else{
        return b;
    }
}

double max_of(double a, double b, double c)
{
    return max_of(max_of(a,b),c);
}

double minmod(double a, double b)
{
    if(a > 0 && b > 0)
    {
        return min_of(a,b);
    }
    else if(a <0 && b < 0)
    {
        return max_of(a,b);
    }
    else{
        return 0;
    }
}

double minmod(double a, double b, double c){
    if(a > 0 && b > 0 && c > 0)
    {
        return min_of(a,b, c);
    }
    else if(a <0 && b < 0 && c < 0)
    {
        return max_of(a,b,c);
    }
    else{
        return 0;
    }
}

double maxmod(double a, double b)
{
    if(a > 0 && b > 0)
    {
        return max_of(a,b);
    }
    else if(a <0 && b < 0)
    {
        return min_of(a,b);
    }
    else{
        return 0;
    }
}

double maxmod(double a, double b, double c){
    if(a > 0 && b > 0 && c > 0)
    {
        return max_of(a,b, c);
    }
    else if(a <0 && b < 0 && c < 0)
    {
        return min_of(a,b,c);
    }
    else{
        return 0;
    }
}
