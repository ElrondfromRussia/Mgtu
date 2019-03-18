#pragma once
#include <omp.h>
#include <cstdlib>
#include <math.h>
#include <ctime>
#include <fstream>
#include "mpi.h"
#include <Windows.h>

#include "matrix.h"

int count_displs(int *sizes, int ran);
Matrix metSoprGr(int size, int rank, Matrix A, Matrix B, int N);
Matrix metZeid(int size, int rank, Matrix A, Matrix B, int N);