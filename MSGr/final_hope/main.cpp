#include "global.h"
//#include "matrix_csr.h"
//#include "mvector.h"


int main(int argc, char** argv)
{
	int rank,size;
	int N = 3;

	double start_time, fin_time; // TIME

	std::vector<double> padiag;
	std::vector<double> paelem;
	//double *padiag = new double[N]; //
	//double *paelem = new double[3]; //
	std::vector<int> paj;
	std::vector<int> pai;
	//int *paj = new int[3];
	//int *pai = new int[N + 1];
	padiag.push_back(4);
	padiag.push_back(2);
	padiag.push_back(4);

	paelem.push_back(1);
	paelem.push_back(2);
	paelem.push_back(4);

	paj.push_back(1);
	paj.push_back(0);
	paj.push_back(1);

	pai.push_back(0);
	pai.push_back(1);
	pai.push_back(1);
	pai.push_back(4);

	std::vector<double> result;
	//double *result = new double[N];
	result.push_back(2);
	result.push_back(1);
	result.push_back(3);


	MPI_Init(&argc, &argv); // START PARALL
	MPI_Comm_size(MPI_COMM_WORLD, &size);	
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	MPI_Barrier(MPI_COMM_WORLD);
	start_time = MPI_Wtime(); // TIME OF START

	Vector Vresult;
	Vresult = sopr_grad_meth(size, rank, paelem, padiag, pai, paj, N, result);

	MPI_Barrier(MPI_COMM_WORLD);
	fin_time = MPI_Wtime() - start_time;
	if (rank == 0)
	{
		std::cout << "\nRESULT TIME = " << fin_time << "\n";
		std::cout << "\nRESULT= " << Vresult << "\n";
	}
	MPI_Finalize(); // FINISH PARALL


	return 0;
}