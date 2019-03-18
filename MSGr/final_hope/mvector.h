#pragma once
#include "global.h"

class Vector
{
	//double *content; //элементы
	std::vector<double> content;
	int N; //размерность

public:
	Vector(); //Конструктор по умолчанию 1
	Vector(int pN); //Конструктор
	Vector(std::vector<double> pcontent, int pN); //конструктор с параметрами
	Vector(const Vector & a); //копирующий конструктор
	~Vector(); //деструктор

	void set_el(int i, double val);
	int get_N() { return N; }

	Vector & operator =(const Vector & b); //перегрузка присваивания
	double operator[](int ind); //перегрузка индексирования
	Vector operator +(Vector & b); //перегрузка сложения
	Vector operator -(Vector & b); //перегрузка вычитания
	Vector operator -(); //перегрузка унарного минуса(результат = -А)
	double operator *(Vector & b); //перегрузка умножения
	Vector operator *(const double & c); //перегрузка умножения

	friend std::ostream & operator <<(std::ostream & a, const Vector & p_a); // вывод содержимого на экран
};
