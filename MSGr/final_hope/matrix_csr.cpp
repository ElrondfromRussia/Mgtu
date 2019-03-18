#include "matrix_csr.h"

///////////////////////////////////////////////////////////////////////////////////////
Matrix_csr::Matrix_csr() // онструктор по умолчанию
{ 
	N = 1;
}

Matrix_csr::Matrix_csr(int pN)
{
	N = pN;	
	for (int i = 0; i < N + 1; i++)
		aI.push_back(0);
	for (int i = 0; i < N ; i++)
		adiag.push_back(0);
}


Matrix_csr::Matrix_csr(std::vector<double> paelem, std::vector<double> padiag, std::vector<int> pai, std::vector<int> paj, int pN, int pN_elem) // constr with parametrs
{
	N = pN;
	N_elem = pN_elem;

	adiag.insert(aelem.begin(), padiag.begin(), padiag.end());
	aelem.insert(aelem.begin(), paelem.begin(), paelem.end());
	aJ.insert(aJ.begin(), paj.begin(), paj.end());
	aI.insert(aI.begin(), pai.begin(), pai.end());

	//std::cout << "\nELEM and J\n";
	//for (int i = 0; i < aelem.size(); i++)
	//{
	//	std::cout << aelem[i] << "  ";
	//}
	//std::cout << "\n";
	//for (int i = 0; i < aJ.size(); i++)
	//{
	//	std::cout << aJ[i] << "  ";
	//}
	//std::cout << "\nDIAG and I\n";
	//for (int i = 0; i < adiag.size(); i++)
	//{
	//	std::cout << adiag[i] << "  ";
	//}
	//std::cout << "\n";
	//for (int i = 0; i < aI.size(); i++)
	//{
	//	std::cout << aI[i] << "  ";
	//}
}

Matrix_csr::Matrix_csr(const Matrix_csr & a) // copy constr
{
	N = a.N;
	N_elem = a.N_elem;

	aelem.insert(aelem.begin(), a.aelem.begin(), a.aelem.end());
	adiag.insert(adiag.begin(), a.adiag.begin(), a.adiag.end());
	aJ.insert(aJ.begin(), a.aJ.begin(), a.aJ.end());
	aI.insert(aI.begin(), a.aI.begin(), a.aI.end());
}

Matrix_csr::~Matrix_csr() //деструктор
{

}
///////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////
Matrix_csr & Matrix_csr::operator =(const Matrix_csr & b)//перегрузка присваивани€
{
	if (this == &b) { return (*this); }
	N = b.N;
	N_elem = b.N_elem;

	aelem.clear();
	adiag.clear();
	aJ.clear();
	aI.clear();

	aelem.insert(aelem.begin(), b.aelem.begin(), b.aelem.end());
	adiag.insert(adiag.begin(), b.adiag.begin(), b.adiag.end());
	aJ.insert(aJ.begin(), b.aJ.begin(), b.aJ.end());
	aI.insert(aI.begin(), b.aI.begin(), b.aI.end());
	
	return(*this);
}

std::ostream & operator <<(std::ostream & p_Out, Matrix_csr & p_M) // вывод содержимого на экран
{
	for (int i = 0; i < p_M.getsizeN(); i++)
	{
		for (int j = 0; j < p_M.getsizeN(); j++)
		{
			printf("%f  ", p_M.get_element(i, j));
		}
		printf("\n");
	}
	printf("\n");
	return p_Out;
}


///////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////
Matrix_csr Matrix_csr::operator +(Matrix_csr & b) //перегрузка сложени€ матриц   $$$$$$$$$$$$$$$$$
{
	omp_set_num_threads(N_TH);

	double k_buf;
	
	std::vector<double> rdiags;
	std::vector<double> raelem;
	std::vector<int> raj;
	std::vector<int> rai;

	int i, j;
	int counter = 0;
	if (b.N == N)
	{
		
#pragma omp parallel private(i,j) shared(b)
		{
#pragma omp for schedule(guided,1)
			for (i = 0; i < N; i++)
			{
				rai.push_back(counter);
				//counter = 0;
				for (j = i; j < N; j++)
				{
					//реализовать заполнение всех массивов TO
					k_buf = get_element(i, j) + b.get_element(i, j);
					if (i == j)
					{
						rdiags.push_back(k_buf);
					}
					if (k_buf != 0)
					{
						if (i != j)
						{
							raelem.push_back(k_buf);
							counter++;
							raj.push_back(j);
						}
					}
				}

			}
		}
		
		Matrix_csr TO(raelem, rdiags, rai, raj, rdiags.size(), raelem.size());
		return TO;
	}
	else
	{
		std::cout << "\n pls Sizes are not suitable\n";
		return Matrix_csr();
	}
}

Matrix_csr Matrix_csr::operator -(Matrix_csr & b) //перегрузка вычитани€ матриц $$$$$$$$$$$$$$$$$
{
	omp_set_num_threads(N_TH);

	double k_buf;
	std::vector<double> rdiags;
	std::vector<double> raelem;
	std::vector<int> raj;
	std::vector<int> rai;

	int i, j;
	int counter = 0;
	if (N == b.N)
	{	
		for ( i = 0; i < 1; i++)
		{
			for (j = 0; j < 1; j++)
			{
				rai.push_back(counter);
				counter = 0;
				//реализовать заполнение всех массивов TO
				k_buf = get_element(i, j) - b.get_element(i, j);
				if (k_buf != 0)
				{
					if (i == j)
					{
						rdiags.push_back(k_buf);
					}
					else
					{
						raelem.push_back(k_buf);
						counter++;
						raj.push_back(j);
					}
				}
			}
		}

		Matrix_csr TO(raelem, rdiags, rai, raj, rdiags.size(), raelem.size());
		return TO;
	}
	else
	{
		std::cout << "\n mns Sizes are not suitable\n";
		return Matrix_csr();
	}
}

Matrix_csr Matrix_csr::operator -() //перегрузка унарного минуса(результат = -ј)
{
	Matrix_csr result(*this);
	int i = 0;
#pragma omp parallel private(i) shared(result)
	{
#pragma omp for schedule(guided,1)
		for (i = 0; i < N; i++)
			result.adiag[i] = -adiag[i];
		for (i = 0; i < N_elem; i++)
			result.aelem[i] = -aelem[i];
	}
	return result;
}

Matrix_csr Matrix_csr::operator *(const double & c) //перегрузка умножени€  $$$$$$$$$$$$$$$$$
{
	omp_set_num_threads(N_TH);
	int i;
	Matrix_csr TO(*this);
#pragma omp parallel private(i) shared(c, TO)
	{
#pragma omp for schedule(guided,1)
		for (i = 0; i < N; i++)
		{
				TO.adiag[i] = TO.adiag[i]*c;
		}
		for (i = 0; i < N_elem; i++)
		{
			TO.aelem[i] = TO.aelem[i] * c;
		}
	}
	return TO;
}
//////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////
Vector Matrix_csr::operator *(Vector & p_b)
{
	omp_set_num_threads(N_TH);
	int i, j;

	if (p_b.get_N() != N)
	{
		std::cout << "\n scal Sizes are not suitable\n";
		std::cout.flush();
		return Vector();
	}
	else
	{		
		//реализовать умножение с заполнением всех массивов
		Vector scal_m(N);
		double buf_sum = 0;
#pragma omp parallel private(j, buf_sum) shared(scal_m, p_b)
		{
#pragma omp for schedule(guided,1)
			for (i = 0; i < 1; i++)
			{
				buf_sum = 0;
				for (j = 0; j < 1; j++)
					buf_sum += get_element(i, j) * p_b[j];					
				scal_m.set_el(i, buf_sum);
			}
		}

		return scal_m;
	}
	return Vector();
}
//////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////
double Matrix_csr::get_element(int i, int j)
{
	double AA = 0;
	int N1, N2;
	if (i < adiag.size() && j < aelem.size())
	{
		if (i == j)
		{
			AA = adiag[i];
		}
		else
		{
			N1 = aI[i];
			N2 = aI[i + 1];
			for (int k = N1; k < N2; k++)
			{
				if (k < aJ.size())
					if (aJ[k] == j)
					{
						AA = aelem[k];
						break;
					}
			}
		}
	}

	return AA;
}
void Matrix_csr::set_element(int i, int j, double val)
{
	//а черт его знает, как элемент добавить
	if (i == j)
	{
		//if (adiag.size() <= i)
		//	adiag.push_back(val);
		//else
		(*this).adiag[i] = val;
	}
	else
	{
		if (val == 0) // val == 0, тогда:
		{
			if (get_element(i, j) != 0 && aI.size() > i) // удалить элемент из aelem (и ai, aj)
			{
				int N1, N2;
				double AA = 0;
				N1 = aI[i];
				N2 = aI[i + 1];
				for (int k = N1; k < N2; k++)
				{
					if (k < N_elem)
						if (aJ[k] == j)
						{
							aJ.erase(aJ.begin() + k);
							aelem.erase(aelem.begin() + k);
							aI[i + 1]--;
							break;
						}
				}
			}
			else  
			{
				// ничего не делать, если такого элемента нет в массивах
			}
		}
		else // val != 0, тогда:
		{
			if (get_element(i, j) != 0 && aI.size() > i) // модифицировать, если такой элемент есть в массивах	
			{
				int N1, N2;
				double AA = 0;
				N1 = aI[i];
				N2 = aI[i + 1];
				for (int k = N1; k < N2; k++)
				{
					if (k < N_elem)
						if (aJ[k] == j)
						{
							aelem[k] = val;
							break;
						}
				}
			}
			else  // добавить элемент, если такого нет в массивах
			{
				int N1, N2;
				double AA = 0;
				N1 = aI[i];
				N2 = aI[i + 1];
				for (int k = N1; k < N2; k++)
				{
					if (k < N_elem)
						if (aJ[k] == j)
						{
							aJ.insert(aJ.begin() + k, j);
							aelem.insert(aelem.begin() + k, val);
							aI[i + 1]++;
							break;
						}
				}
			}
		}
	}
	std::cout << *this;
}
//////////////////////////////////////////////////////////////////////////////////////

int Matrix_csr::getsizeN()
{
	return N;
}

int Matrix_csr::getsizeN_elem()
{
	return N_elem;
}

/////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////