#pragma once
#include "global.h"

class Vector
{
	//double *content; //��������
	std::vector<double> content;
	int N; //�����������

public:
	Vector(); //����������� �� ��������� 1
	Vector(int pN); //�����������
	Vector(std::vector<double> pcontent, int pN); //����������� � �����������
	Vector(const Vector & a); //���������� �����������
	~Vector(); //����������

	void set_el(int i, double val);
	int get_N() { return N; }

	Vector & operator =(const Vector & b); //���������� ������������
	double operator[](int ind); //���������� ��������������
	Vector operator +(Vector & b); //���������� ��������
	Vector operator -(Vector & b); //���������� ���������
	Vector operator -(); //���������� �������� ������(��������� = -�)
	double operator *(Vector & b); //���������� ���������
	Vector operator *(const double & c); //���������� ���������

	friend std::ostream & operator <<(std::ostream & a, const Vector & p_a); // ����� ����������� �� �����
};
