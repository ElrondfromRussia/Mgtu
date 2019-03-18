#include "matrix.h"
#include <ctime>
#include <cstdlib>
#include <omp.h>

#define N_TH 1 // количество процессов для OpenMP

///////////////////////////////////////////////////////////////////////////////////////
Matrix::Matrix() //Конструктор по умолчанию
{ 
	Cont = NULL;
	n = 0;
	m = 0;
}


Matrix::Matrix(int N, int M) // constr with parametrs
{
	n = N;
	m = M;
	Cont = new double *[n];
	for (int i = 0; i < n; i++)
	{
		Cont[i] = new double[m];
	}
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
		{
			Cont[i][j] = 0;
		}
	}
}

Matrix::Matrix(const Matrix & a) // copy constr
{
	n = a.n;
	m = a.m;
	Cont = new double *[n];
	{
		for (int i = 0; i < n; i++)
		{
			Cont[i] = new double[m];
		}
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < m; j++)
			{
				Cont[i][j] = a.Cont[i][j];
			}
		}
	}
}

Matrix::~Matrix() //деструктор
{
	if (Cont)
	{
		for (int i = 0; i < n; i++)
			delete[] Cont[i];
	}
	delete[]Cont;
}
///////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////
int Matrix::fill() //  fill
{
	//позже тут необходимо реализовать заполнение матриц из файла

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
		{
			
		}
	}
	return 0;
}
///////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////
void Matrix::rand_fill(const int p_num) // random fill
{
	srand(time(NULL));

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
		{
			Cont[i][j] = (1 + rand() % 10) / 2.0 + p_num*(1 + rand() % 10);
		}
	}
}
///////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////
Matrix & Matrix::operator =(const Matrix & b)//перегрузка присваивания
{
	if (this == &b) { return (*this); }
	if (Cont)
	{	
		for (int j = 0; j < n; j++)
			delete[] Cont[j];
	}
	delete []Cont;
	n = b.n;
	m = b.m;
	Cont = new double *[n];
	for (int i = 0; i < n; i++)
	{
		Cont[i] = new double[m];
	}
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
		{
			Cont[i][j] = b.Cont[i][j];
		}
	}
	return(*this);
}

std::ostream & operator <<(std::ostream & p_Out, const Matrix & p_M) // вывод содержимого на экран
{
	for (int i = 0; i < p_M.n; i++)
	{
		for (int j = 0; j < p_M.m; j++)
			printf("%f  ", p_M.Cont[i][j]);
		printf("\n");
	}
	printf("\n");
	return p_Out;
}


///////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////
Matrix Matrix::transp() // транспонирование матрицы
{
	Matrix c(m, n);
	for (int i = 0; i < n; i++)
	for (int j = 0; j < m; j++)
	{
		c.Cont[j][i] = Cont[i][j];
	}
	return c;
}

int Matrix::sim_matr() // simmetrir                    $$$$$$$$$$$$$$$$$
{
	int i, j;
	omp_set_num_threads(N_TH);
	printf("\n we have %d OpenMP threads\n", omp_get_num_threads());
	if (n == m)
	{
#pragma omp parallel private(i,j)
		{
#pragma omp for schedule(guided,1)
			for (i = 0; i < n; i++) // симметрировали
			{
				//printf("\n %d thread\n", omp_get_thread_num());
				for (j = i; j < m; j++)
				{
					Cont[j][i] = Cont[i][j];
					//if (i == j)
					//	Cont[j][i] *= 10;
				}
			}
		}		
	}
	else
	{
		std::cout << "\n Not squere matrix\n";
		return -1;
	}
	return 0;
}
///////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////
double * Matrix::operator [](int Str) //перегрузка индексирования
{
	return Cont[Str];
}
///////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////
// Получение матрицы без i-й строки и j-го столбца
void GetMatr(Matrix & mas, Matrix & p, int i, int j, int m) {
	int ki, kj, di, dj;
	di = 0;
	for (ki = 0; ki<m - 1; ki++) { // проверка индекса строки
		if (ki == i) di = 1;
		dj = 0;
		for (kj = 0; kj<m - 1; kj++) { // проверка индекса столбца
			if (kj == j) dj = 1;
			p[ki][kj] = mas[ki + di][kj + dj];
		}
	}
}

// Рекурсивное вычисление определителя
double Determinant(Matrix & mas, int m) {
	int i, j, n;
	double d, k;
	Matrix p(m, m);

	j = 0; d = 0;
	k = 1;
	n = m - 1;
	if (m<1) std::cout << "Can't calculate determinant!";
	if (m == 1) {
		d = mas[0][0];
		return(d);
	}
	if (m == 2) {
		d = mas[0][0] * mas[1][1] - (mas[1][0] * mas[0][1]);
		return(d);
	}
	if (m>2) {
		for (i = 0; i<m; i++) {
			GetMatr(mas, p, i, 0, m);
			d = d + k * mas[i][0] * Determinant(p, n);
			k = -k;
		}
	}
	return(d);
}
//////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////
Matrix Matrix::operator +(const Matrix & b) //перегрузка сложения матриц   $$$$$$$$$$$$$$$$$
{
	omp_set_num_threads(N_TH);

	int i, j;
	if ((n == b.n)&(m == b.m))
	{
		Matrix TO(n, m);
#pragma omp parallel private(i,j) shared(TO,b)
		{
#pragma omp for schedule(guided,1)
			for (i = 0; i < n; i++)
			{
				for (j = 0; j < m; j++)
				{
					TO.Cont[i][j] = Cont[i][j] + b.Cont[i][j];
				}
			}
		}
		return TO;
	}
	else
	{
		std::cout << "\n pls Sizes are not suitable\n";
		Matrix TO(0, 0);
		return TO;
	}
}

Matrix Matrix::operator -(const Matrix & b) //перегрузка вычитания матриц $$$$$$$$$$$$$$$$$
{
	int i, j;
	omp_set_num_threads(N_TH);
	if ((n == b.n)&(m == b.m))
	{
		Matrix TO(n, m);
#pragma omp parallel private(i,j) shared(TO,b)
		{
#pragma omp for schedule(guided,1)
			for ( i = 0; i < n; i++)
			{
				for (j = 0; j < m; j++)
				{
					TO.Cont[i][j] = Cont[i][j] - b.Cont[i][j];
				}
			}
		}
		return TO;
	}
	else
	{
		std::cout << "\n mns Sizes are not suitable\n";
		Matrix TO(0, 0);
		return TO;
	}
}

Matrix Matrix::operator *(const Matrix & b) //перегрузка умножения  $$$$$$$$$$$$$$$$$
{
	omp_set_num_threads(N_TH);
	int i, j, k;
	if (m == b.n)
	{
		Matrix result(n, b.m);
#pragma omp parallel private(j,k) shared(result,b)
		{
#pragma omp for schedule(guided,1)
			for (i = 0; i < n; i++)
			{
				for (j = 0; j < b.m; j++)
				{
					result.Cont[i][j] = 0;
					for (k = 0; k < b.n; k++)
						result.Cont[i][j] += this->Cont[i][k] * b.Cont[k][j];
				}
			}
		}
		return result;
	}
	else
	{
		std::cout << "\nMULT_MATR size error\n";
		Matrix result(0, 0);
		return result;
	}
}

Matrix Matrix::operator -() //перегрузка унарного минуса(результат = -А)
{
	/*Matrix result(n, m);
	for (int i = 0; i < n; i++)
	for (int j = 0; j < m; j++)
	{
		result.Cont[i][j] = -Cont[i][j];
	}
	return result;*/
	for (int i = 0; i < n; i++)
		for (int j = 0; j < m; j++)
		{
			Cont[i][j] = -Cont[i][j];
		}
	return (*this);
}

Matrix Matrix::operator *(const double & c) //перегрузка умножения  $$$$$$$$$$$$$$$$$
{
	omp_set_num_threads(N_TH);
	int i, j;
	Matrix TO(n, m);
#pragma omp parallel private(i,j) shared(c, TO)
	{
#pragma omp for schedule(guided,1)
		for (i = 0; i < n; i++)
		{
			for (j = 0; j < m; j++)
			{
				TO.Cont[i][j] = c*Cont[i][j];
			}
		}
	}
	return TO;
}
//////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////
double  scal_mult(Matrix & p_a, Matrix & p_b)
{
	if ((p_a.n != 1) || (p_b.m != 1) || (p_a.m != p_b.n))
	{
		std::cout << "\n scal Sizes are not suitable\n";
		return 1;
	}
	else
	{
		Matrix res = p_a*p_b;
		double scal_m = res[0][0];
		return scal_m;
	}
}
//////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////
double get_element(int i, int j, double *adiag, double *altr, int* jptr, int* iptr)
{
	double AA = 0;
	int N1, N2;

	if (i >= j)
	{
		N1 = iptr[i - 1];
		N2 = iptr[i];
		if (i == j)
			AA = adiag[i - 1];
		else
		{
			for (int k = N1; k < N2; k++)
			{
				if (jptr[k - 1] == j)
				{
					AA = altr[k - 1];
					return AA;
				}
			}
		}
	}
	else if (i<j)
	{
		N1 = iptr[j - 1];
		N2 = iptr[j];
		for (int k = N1; k < N2; k++)
		{
			if (jptr[k - 1] == i)
			{
				AA = altr[k - 1];
				return AA;
			}
		}
	}
	return AA;
}
//////////////////////////////////////////////////////////////////////////////////////