#pragma once

#include "global.h"
#include "mvector.h"


class Matrix_csr
{
	//double *aelem; //�������� �� �� ���������
	std::vector<double> aelem;
	//double *adiag; // ������������ ��������
	std::vector<double> adiag;
	//int *aI; //������ 
	std::vector<int> aI;
	//int *aJ; //�������
	std::vector<int> aJ;
	int N;  //���������� ���������
	int N_elem;

public:
	Matrix_csr(); //����������� �� ���������
	Matrix_csr(int pN); //����������� 
	Matrix_csr(std::vector<double> paelem, std::vector<double> padiag, std::vector<int> pai, std::vector<int> paj, int pN, int pN_elem); //����������� � �����������
	Matrix_csr(const Matrix_csr & a); //���������� �����������
	~Matrix_csr(); //����������

	Matrix_csr & operator =(const Matrix_csr & b); //���������� ������������

	friend std::ostream & operator <<(std::ostream & a, Matrix_csr & p_a); // ����� ����������� �� �����

	double get_element(int i, int j); //���������� ��������������
	void set_element(int i, int j, double val); //���������� ��������������
	int getsizeN();
	int getsizeN_elem();

	Matrix_csr operator +(Matrix_csr & b); //���������� �������� ������
	Matrix_csr operator -(Matrix_csr & b); //���������� ��������� ������
	Matrix_csr operator -(); //���������� �������� ������(��������� = -�)
	Matrix_csr operator *(const double & c); //���������� ���������
	Vector operator *(Vector & p_b); //���������� ���������
};


