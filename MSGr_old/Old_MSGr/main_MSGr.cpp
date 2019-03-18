#include "global.h"


//mpiexec -np 1 Old_MSGr.exe

const int N = 3;

int main(int argc, char** argv)
{
	int rank, size; // for parallel realisation

	double start_time, fin_time; // TIME

	MPI_Init(&argc, &argv); // START PARALL
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	MPI_Barrier(MPI_COMM_WORLD);
	start_time = MPI_Wtime(); // TIME OF START

	////////////////////////////////////////////////////

	/////////////////тестовые для А
	double *padiag = new double[N]; //
	double *paelem = new double[3]; //
	int *paj = new int[3];
	int *pai = new int[N + 1];	
	padiag[0] = 4;
	padiag[1] = 2;
	padiag[2] = 4;

	paelem[0] = 1;
	paelem[1] = 2;
	paelem[2] = 4;

	paj[0] = 1;
	paj[1] = 0;
	paj[2] = 1;

	pai[0] = 0;
	pai[1] = 1;
	pai[2] = 1;
	pai[3] = 4;

	Matrix A(N, N); // Матрица СЛАУ
	A.read_from_massives(padiag, paelem, paj, pai); //читаем матрицу СЛАУ

	//////////////////тестовые для B
	double *bpmass = new double[N]; //

	bpmass[0] = 2;
	bpmass[1] = 1;
	bpmass[2] = 3;

	Matrix B(N, 1); //правый столбец
	B.read_column(bpmass); //читаем столбец правый

	//сразу чистим память из-под массивов
	delete[]padiag;
	delete[]paelem;
	delete[]paj;
	delete[]pai;

	//печатаем тестовое слау
	printf("\nrank = %d\n", rank);
	if (rank == 0)
	{
		std::cout << "\na:\n" << A;
		std::cout << "\nb:\n" << B;
	}
	////////////////////////////////////////////////////

	Matrix ReS(N, 1);//Сюда принимаем решение

	///////////САМ АЛГОРИТМ
	ReS = metSoprGr(size, rank, A, B, N);
	//////////КОНЕЦ (позже будут другие варианты)


	MPI_Barrier(MPI_COMM_WORLD); // END_TIME
	fin_time = MPI_Wtime() - start_time;
	if (rank == 0)
	{
		std::cout << "\nRESULT TIME = " << fin_time << "\n";
		std::cout << "\nSOLUTION:\n";
		std::cout << ReS;
	}

	MPI_Finalize(); // FINISH PARALL
	return 0;
}