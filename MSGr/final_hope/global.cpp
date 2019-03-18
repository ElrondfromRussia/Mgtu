#include "global.h"

int count_displs(int *sizes, int ran)
{
	int disp = 0;
	for (int h = 0; h < ran; h++)
		disp += sizes[h];
	return disp;
}

Vector sopr_grad_meth(int size, int rank, std::vector<double> paelem, std::vector<double> padiag, std::vector<int> pai, std::vector<int> paj, int pN, std::vector<double> pb_col)
{
	int i, j;
	int n, jbeg, jend; // for parallel realisation

					   ///////////////////////////////////////////////////////////////////////////
	Matrix_csr A(paelem, padiag, pai, paj, pN, sizeof(paelem)); // ������� ����

	Vector B(pb_col, pN); // ������� ��������� ������. ����� ������� ��� ���� ����-�

	printf("\nrank = %d\n", rank);
	if (rank == 0)
	{
		std::cout << "\na:\n" << A;
		std::cout << "\nb:\n" << B;
	}

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// ��� ������������ ���������/������ ������� � ������� ��� ���������� �� ������ �����
	n = (pN - 1) / size + 1; // for deviding matrix A
	jbeg = rank*n;
	jend = (rank + 1)*n;

	int *sizes = new int[size];
	int *displs = new int[size];

	int s_b = ((jend > pN) ? pN : jend) - jbeg;

	MPI_Gather(&s_b, 1, MPI_INT, sizes, 1, MPI_INT, 0, MPI_COMM_WORLD);

	if (rank == 0)
	{
		for (int k = 0; k < size; k++)
		{
			displs[k] = count_displs(sizes, k);
		}
	}

	n = s_b; // ������ n �� �����, ��� �������� ������

	Matrix_csr A_part(n); // ������� ���� partial
	Vector B_part(n); // ������� ��������� ������ ���������

	for (j = jbeg; j < ((jend > pN) ? pN : jend); j++)
	{
		B_part.set_el(j%n, B[j]);
		for (int k = 0; k < pN; k++)
		{
			std::cout << "\ni: " << j%n << "   j: " << k << "  val:" << A.get_element(j, k) << " \n";
			std::cout.flush();
			A_part.set_element(j%n, k, A.get_element(j, k));
			std::cout << "\nApartel: \n" << A_part.get_element(j%n, k);
			std::cout.flush();
		}
	}
	
	printf("\nrank = %d\n", rank);
	if (rank == 0)
	{
		std::cout << "\napart:\n" << A_part;
		std::cout.flush();
		std::cout << "\nbpart:\n" << B_part;
		std::cout.flush();
	}
	
	Vector x1(pN);
	Vector x2(pN); // � �� ����� �� �����, ��� ��� ������ ���� ������������ ��� ��������� �� ��������� ������� ����

	Vector d1(pN);
	Vector d2(pN); // ������� d ���� �� �����, ��� ����������� ��������� ����� ����� 

	Vector pre_alf1(pN);
	Vector pre_alf2(pN); // pre_alf

	Vector g_part1 = B_part;
	Vector g_part2(n);  // ��������� ������� g

	Vector g1 = B;
	Vector g2(pN);  //  ������� g

	Vector part_pre_alf1(n);
	Vector part_pre_alf2(n); // ��������� ������� g

	double *alpha = new double[pN + 1];
	
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// ��� ����� ���������������� ��� �������� 
	// � ������ ����� ����������� �������, ������� ����� ���������� � 0-� ����� � ������ ��������

	double *small_buf = new double[n];
	double *big_buf = new double[pN];

	for (i = 0; i < pN; i++) ///////////////////////////////////////
	{
		MPI_Barrier(MPI_COMM_WORLD);
		// for g
		g_part2 = A_part*x1 - B_part;
		
		for (int y = 0; y < n; y++)
		{
			small_buf[y] = g_part2[y];
		}

		MPI_Gatherv(small_buf, s_b, MPI_DOUBLE, big_buf, sizes, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

		

		//
		// for d
		if (rank == 0)
		{
			for (int y = 0; y < pN; y++)
			{
				g2.set_el(y, big_buf[y]);
			}
			d2 = -g2 + d1 * ((g2 * g2) / (g2 * g1));
			for (int y = 0; y < pN; y++)
				big_buf[y] = d2[y];
		}

		MPI_Bcast(big_buf, pN, MPI_DOUBLE, 0, MPI_COMM_WORLD);

		if (rank != 0)
		{
			for (int y = 0; y < pN; y++)
				d2.set_el(y, big_buf[y]);
		}
		//
		//for alpha
		part_pre_alf2 = A_part*d2;
		for (int y = 0; y < n; y++)
			small_buf[y] = part_pre_alf2[y];

		MPI_Gatherv(small_buf, s_b, MPI_DOUBLE, big_buf, sizes, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

		if (rank == 0)
		{
			for (int y = 0; y < pN; y++)
				pre_alf2.set_el(y, big_buf[y]);
			alpha[i + 1] = (d2 * g2) / (d2 * pre_alf2);
			//
			// for x
			x2 = x1 + d2 * alpha[i + 1];
			for (int y = 0; y < pN; y++)
				big_buf[y] = x2[y];
		}

		MPI_Bcast(big_buf, pN, MPI_DOUBLE, 0, MPI_COMM_WORLD);

		if (rank != 0)
		{
			for (int y = 0; y < pN; y++)
				x2.set_el(y, big_buf[y]);
		}
		//
		//std::cout << "\nITER END: " << i;			
		x1 = x2;
		d1 = d2;
		g1 = g2;
	}
	
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// ��� ������ �� ��� ���� ������� �������
	delete[]small_buf;
	delete[]big_buf;
	delete[]alpha;
	delete[]displs;
	delete[]sizes;
	///////////////////////////////////////////////////////////////////////////////
	return 0;
}