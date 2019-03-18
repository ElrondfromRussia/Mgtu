#include "global.h"


int count_displs(int *sizes, int ran)
{
	int disp = 0;
	for (int h = 0; h < ran; h++)
		disp += sizes[h];
	return disp;
}

Matrix metSoprGr(int size, int rank, Matrix A, Matrix B, int N)
{
	int i, j;
	int n, jbeg, jend; // for parallel realisation

	n = (N - 1) / size + 1; // for deviding matrix 
	jbeg = rank*n;
	jend = (rank + 1)*n;

	int *sizes = new int[size];
	int *displs = new int[size];

	int s_b = ((jend > N) ? N : jend) - jbeg;

	std::cout << "\nrank: " << rank << "    s_b: " << s_b;
	std::cout.flush();

	MPI_Gather(&s_b, 1, MPI_INT, sizes, 1, MPI_INT, 0, MPI_COMM_WORLD);

	if (rank == 0)
	{
		for (int k = 0; k < size; k++)
		{
			displs[k] = count_displs(sizes, k);
		}
	}

	n = s_b; // меняем n на новое, для неравных кусков

	Matrix A_part(n, N); // Матрица СЛАУ partial
	Matrix B_part(n, 1); // Столбец свободных членов частичный

	for (j = jbeg; j < ((jend > N) ? N : jend); j++)
	{
		B_part[j%n][0] = B[j][0];
		for (int k = 0; k < N; k++)
			A_part[j%n][k] = A[j][k];
	}

	Matrix x1(N, 1);
	Matrix x2(N, 1); // х Не делим на куски, так как должно быть соответствие при умножении на частичную матрицу СЛАУ

	Matrix d1(N, 1);
	Matrix d2(N, 1); // столбец d тоже не делим, для размерности умножения строк нужно 

	Matrix pre_alf1(N, 1);
	Matrix pre_alf2(N, 1); // pre_alf

	Matrix g_part1 = B_part;
	Matrix g_part2(n, 1);  // Частичный столбец g

	Matrix g1 = B;
	Matrix g2(N, 1);  //  столбец g

	Matrix part_pre_alf1(n, 1);
	Matrix part_pre_alf2(n, 1); // Частичный столбец g

	double *alpha = new double[N + 1];

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Тут будем распараллеливать сам алгоритм 
	// В каждом проце вычисляются кусочки, которые потом собираются в 0-м проце в каждой итерации

	double *small_buf = new double[n];
	double *big_buf = new double[N];

	for (i = 0; i < N; i++) ///////////////////////////////////////
	{
		MPI_Barrier(MPI_COMM_WORLD);
		// for g
		g_part2 = A_part*x1 - B_part;
		for (int y = 0; y < n; y++)
		{
			small_buf[y] = g_part2[y][0];
		}

		MPI_Gatherv(small_buf, s_b, MPI_DOUBLE, big_buf, sizes, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

		//
		// for d
		if (rank == 0)
		{
			for (int y = 0; y < N; y++)
			{
				g2[y][0] = big_buf[y];
			}
			d2 = -g2 + d1 * ((scal_mult(g2.transp(), g2)) / (scal_mult(g1.transp(), g1)));
			for (int y = 0; y < N; y++)
				big_buf[y] = d2[y][0];
		}

		MPI_Bcast(big_buf, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

		if (rank != 0)
		{
			for (int y = 0; y < N; y++)
				d2[y][0] = big_buf[y];
		}
		//
		//for alpha
		part_pre_alf2 = A_part*d2;
		for (int y = 0; y < n; y++)
			small_buf[y] = part_pre_alf2[y][0];

		MPI_Gatherv(small_buf, s_b, MPI_DOUBLE, big_buf, sizes, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

		if (rank == 0)
		{
			for (int y = 0; y < N; y++)
				pre_alf2[y][0] = big_buf[y];
			alpha[i + 1] = (scal_mult(d2.transp(), g2)) / (scal_mult(d2.transp(), pre_alf2));
			//
			// for x
			x2 = x1 + d2 * alpha[i + 1];
			for (int y = 0; y < N; y++)
				big_buf[y] = x2[y][0];
		}

		MPI_Bcast(big_buf, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

		if (rank != 0)
		{
			for (int y = 0; y < N; y++)
				x2[y][0] = big_buf[y];
		}
		//
		//std::cout << "\nITER END: " << i;			
		x1 = x2;
		d1 = d2;
		g1 = g2;
	}
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// чистилище
	delete[]small_buf;
	delete[]big_buf;
	delete[]alpha;
	delete[]displs;
	delete[]sizes;

	return x2;
}

Matrix metZeid(int size, int rank, Matrix A, Matrix B, int N)
{
	Matrix res;
	return res;
}