#include "Distributions.h"
#include <random>
#include <cmath>

// Функция для генерации случайного числа в диапазоне [0, 1)
double random() {
	return (double)rand() / (RAND_MAX + 1.0);
}

//--Отрицательное биномиальное распределение 1--
/*
Args:
	m - количество успехов
	p - вероятность успеха
*/
int NegativeBinomialDistribution1(int m, double p, double r) {
	double q = 1 - p;
	p = pow(p, m);

	//// Приводим к 0 < r < 1
	//while (r >= 1) {
	//	r /= 10;
	//}
	r = random();

	int z = 0;

	while (true) {
		r -= p;
		if (r < 0) {break;}

		z++;
		p *= q * (m - 1 + z) / z;
	}

	return z;
}

//--Распределение Паскаля--
int PascalsDistribution(int m, double p, double r) {
	return NegativeBinomialDistribution1(m, p, r) + m;
}
