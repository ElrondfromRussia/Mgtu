#include "Matrix_csr.h"

#define N_TH 1 // количество процессов для OpenMP
/////////////////////////////////////////////////////////////////////////////////////////
#include "mvector.h"

Vector::Vector()
{
	N = 1;
}

Vector::Vector(int pN)
{
	N = pN;
}

Vector::Vector(std::vector<double> pcontent, int pN)
{
	N = pN;

	for (int i = 0; i < pcontent.size(); i++)
	{
		content.push_back(pcontent[i]);
	}
}

Vector::Vector(const Vector & a)
{
	N = a.N;

	for (int i = 0; i < a.content.size(); i++)
	{
		content.push_back(a.content[i]);
	}
}

Vector::~Vector()
{
}

void Vector::set_el(int i, double val)
{
	if (i < content.size())
		content[i] = val;
	else
		content.push_back(val);
}

Vector & Vector::operator=(const Vector & b)
{
	if (this == &b) { return (*this); }
	N = b.N;
	content.clear();
	for (int i = 0; i < N; i++)
	{
		content.push_back(b.content[i]);
	}
	return(*this);
}

double Vector::operator [](int ind) //перегрузка индексирования
{
	return content[ind];
}

Vector Vector::operator+(Vector & b)
{
	omp_set_num_threads(N_TH);

	int i;
	if (N == b.N)
	{
		Vector TO(N);
#pragma omp parallel private(i) shared(TO,b)
		{
#pragma omp for schedule(guided,1)
			for (i = 0; i < N; i++)
			{
				TO.content[i] = content[i] + b.content[i];
			}
		}
		return TO;
	}
	else
	{
		std::cout << "\n pls Sizes are not suitable\n";
	}
	return Vector();
}

Vector Vector::operator-(Vector & b)
{
	omp_set_num_threads(N_TH);

	int i;
	if (N == b.N)
	{
		Vector TO(N);
#pragma omp parallel private(i) shared(TO,b)
		{
#pragma omp for schedule(guided,1)
			for (i = 0; i < N; i++)
			{
				TO.content[i] = content[i] - b.content[i];
			}
		}
		return TO;
	}
	else
	{
		std::cout << "\n substr Sizes are not suitable\n";
	}
	return Vector();
}

Vector Vector::operator-()
{
	return Vector();
}

double Vector::operator*(Vector & b)
{
	omp_set_num_threads(N_TH);

	int i;
	if (N == b.N)
	{
		double TO = 0;
#pragma omp parallel private(i) shared(TO,b)
		{
#pragma omp for schedule(guided,1)
			for (i = 0; i < N; i++)
			{
				TO += content[i] * b.content[i];
			}
		}
		return TO;
	}
	else
	{
		std::cout << "\n mult Sizes are not suitable\n";
	}
	return 0;
}

Vector Vector::operator*(const double & c)
{
	omp_set_num_threads(N_TH);

	int i;
	Vector TO(N);
#pragma omp parallel private(i) shared(TO)
	{
#pragma omp for schedule(guided,1)
		for (i = 0; i < N; i++)
		{
			TO.content[i] = content[i] * c;
		}
	}
	return TO;
}

std::ostream & operator <<(std::ostream & p_Out, const Vector & p_M) // вывод содержимого на экран
{
	for (int i = 0; i < p_M.N; i++)
	{
		printf("%f  ", p_M.content[i]);
	}
	printf("\n");
	return p_Out;
}

