#pragma once

#include <omp.h>
#include <cstdlib>
#include <math.h>
#include <ctime>
#include <fstream>
#include "mpi.h"
#include <Windows.h>
#include <vector>
#include <iostream>
#include <conio.h>

#include "matrix_csr.h"
#include "mvector.h"

#define N_TH 1 // количество процессов для OpenMP

int count_displs(int *sizes, int ran);
Vector sopr_grad_meth(int size, int rank, std::vector<double> paelem, std::vector<double> padiag, std::vector<int> pai, std::vector<int> paj, int pN, std::vector<double> pb_col);
