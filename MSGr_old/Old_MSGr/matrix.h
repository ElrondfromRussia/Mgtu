#pragma once
#include <iostream>
#include <ctime>
#include <cstdlib>
#include <conio.h>
#include <fstream>

class Matrix
{
	double **Cont; //содержание матрицы
	int n, m;  //количество строк, столбцов

public:
	Matrix(); // онструктор по умолчанию
	Matrix(int N, int M); //конструктор с параметрами
	Matrix(const Matrix & a); //копирующий конструктор
	~Matrix(); //деструктор

	void read_from_massives(double* padiag, double* paelem, int* paj, int* pai);
	void read_column(double* pmass);

	Matrix & operator =(const Matrix & b); //перегрузка присваивани€

	int fill();
	void rand_fill(const int p_num); // random fill

	friend std::ostream & operator <<(std::ostream & a, const Matrix & p_a); // вывод содержимого на экран

	Matrix transp(); // транспонирование матрицы
	int sim_matr(); // симметрирование матрицы

	double*operator[](int St); //перегрузка индексировани€

	friend void GetMatr(Matrix & mas, Matrix & p, int i, int j, int m); // ѕолучение матрицы без i-й строки и j-го столбца
	friend double Determinant(Matrix & mas, int m); // ¬ычисление определител€ (рекурсивное)

	Matrix operator +(const Matrix & b); //перегрузка сложени€ матриц
	Matrix operator -(const Matrix & b); //перегрузка вычитани€ матриц
	Matrix operator -(); //перегрузка унарного минуса(результат = -ј)
	Matrix operator *(const Matrix & b); //перегрузка умножени€
	Matrix operator *(const double & c); //перегрузка умножени€

	friend double scal_mult(Matrix & p_a, Matrix & p_b);	
};