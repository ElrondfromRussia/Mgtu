#pragma once
#include <iostream>
#include <ctime>
#include <cstdlib>
#include <conio.h>
#include <fstream>

class Matrix
{
	double **Cont; //���������� �������
	int n, m;  //���������� �����, ��������

public:
	Matrix(); //����������� �� ���������
	Matrix(int N, int M); //����������� � �����������
	Matrix(const Matrix & a); //���������� �����������
	~Matrix(); //����������

	void read_from_massives(double* padiag, double* paelem, int* paj, int* pai);
	void read_column(double* pmass);

	Matrix & operator =(const Matrix & b); //���������� ������������

	int fill();
	void rand_fill(const int p_num); // random fill

	friend std::ostream & operator <<(std::ostream & a, const Matrix & p_a); // ����� ����������� �� �����

	Matrix transp(); // ���������������� �������
	int sim_matr(); // ��������������� �������

	double*operator[](int St); //���������� ��������������

	friend void GetMatr(Matrix & mas, Matrix & p, int i, int j, int m); // ��������� ������� ��� i-� ������ � j-�� �������
	friend double Determinant(Matrix & mas, int m); // ���������� ������������ (�����������)

	Matrix operator +(const Matrix & b); //���������� �������� ������
	Matrix operator -(const Matrix & b); //���������� ��������� ������
	Matrix operator -(); //���������� �������� ������(��������� = -�)
	Matrix operator *(const Matrix & b); //���������� ���������
	Matrix operator *(const double & c); //���������� ���������

	friend double scal_mult(Matrix & p_a, Matrix & p_b);	
};