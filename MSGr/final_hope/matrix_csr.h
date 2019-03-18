#pragma once

#include "global.h"
#include "mvector.h"


class Matrix_csr
{
	//double *aelem; //элементы не на диагонали
	std::vector<double> aelem;
	//double *adiag; // диагональные элементы
	std::vector<double> adiag;
	//int *aI; //строки 
	std::vector<int> aI;
	//int *aJ; //столбцы
	std::vector<int> aJ;
	int N;  //количество элементов
	int N_elem;

public:
	Matrix_csr(); // онструктор по умолчанию
	Matrix_csr(int pN); // онструктор 
	Matrix_csr(std::vector<double> paelem, std::vector<double> padiag, std::vector<int> pai, std::vector<int> paj, int pN, int pN_elem); //конструктор с параметрами
	Matrix_csr(const Matrix_csr & a); //копирующий конструктор
	~Matrix_csr(); //деструктор

	Matrix_csr & operator =(const Matrix_csr & b); //перегрузка присваивани€

	friend std::ostream & operator <<(std::ostream & a, Matrix_csr & p_a); // вывод содержимого на экран

	double get_element(int i, int j); //перегрузка индексировани€
	void set_element(int i, int j, double val); //перегрузка индексировани€
	int getsizeN();
	int getsizeN_elem();

	Matrix_csr operator +(Matrix_csr & b); //перегрузка сложени€ матриц
	Matrix_csr operator -(Matrix_csr & b); //перегрузка вычитани€ матриц
	Matrix_csr operator -(); //перегрузка унарного минуса(результат = -ј)
	Matrix_csr operator *(const double & c); //перегрузка умножени€
	Vector operator *(Vector & p_b); //перегрузка умножени€
};


