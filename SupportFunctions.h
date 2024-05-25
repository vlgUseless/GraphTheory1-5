#pragma once
#include "Distributions.h"
#include <iostream>
#include <string>
#include <vector>
#include <iomanip>
using namespace std;


class OrientedGraph;

int InputVertices();
int InputVertex(OrientedGraph* OrGraph, int numb);
int InputEdges(OrientedGraph* OrGraph);
bool InputMinOrMax();
bool InputWeightPar();
void PrintList2d(const vector<vector<int>> List, string ListName);

template<typename T>
void PrintList(const vector<T> List, string ListName) {
	cout << "---" << ListName << "---" << endl;
	for (int i = 0; i < List.size(); i++) {
		cout << i << ": " << List[i] << endl;
	}
	cout << endl;
}

template<typename T>
void PrintMatrix(const vector<vector<T>>& Matrix, string MatrixName) {
	cout << setiosflags(ios::left);

	cout << "---" << MatrixName << "---" << endl;
	for (int i = 0; i < Matrix.size(); i++) {
		cout << setw(3) << "|";
		for (int j = 0; j < Matrix.size(); j++) {
			cout << setw(5) << Matrix[i][j];
		}
		cout <<"|" << endl;
	}
	cout << endl;
}

template<typename T>
vector<vector<T>> InputMatrix() {
	int vertices = InputVertices();
	// allocating memory
	vector<vector<T>> Matrix(vertices, vector<T>(vertices, 0));
	for (int i = 0; i < Matrix.size(); i++) {
		for (int j = 0; j < Matrix.size(); j++) {
			string temp;
			cin >> temp;
			if (temp == "inf") {
				Matrix[i][j] = INFINITY;
			}
			else {
				Matrix[i][j] = (double)stoi(temp);
			}
		}
	}

	return Matrix;
}

double GaussDet(vector<vector<double>> m);
void InputWeightMatrix(vector<vector<double>>& Matrix, const vector<vector<int>>& AdjMatrix);
void Menu();
void InputProcessing(string Action, OrientedGraph*& OrGraph);
vector<vector<int>> MultiplyMatrix(const vector<vector<int>>& Matrix1, const vector<vector<int>>& Matrix2);
void UnionMatrix(vector<vector<int>>& Matrix1, const vector<vector<int>>& Matrix2);

void TestShimbelsMethod();

vector<vector<int>> AddMatrix(const vector<vector<int>>& Matrix1, const vector<vector<int>>& Matrix2);

// јвтоматическое заполнение матрицы весов в соответствии с распределением ѕаскал€ с параметрами m и p
template<typename T>
void GenerateMatrix(vector<vector<T>>& Matrix, const vector<vector<int>>& AdjMatrix, int Nedge, int m, double p, bool PosParam) {
	// If weights can be negative
	if (PosParam == 0) {
		// Number of each weights
		int NPositive = PascalsDistribution(Nedge / 2, p, rand());
		int NNegative = Nedge - NPositive;

		int c = 0;
		// Go through every row
		for (int i = 0; i < Matrix.size(); i++) {
			// Go through every column
			for (int j = 0; j < Matrix[i].size(); j++) {
				// If vertices are adjacent
				if (AdjMatrix[i][j] == 1) {
					// Alternate positive and negative weights
					if (c % 2 == 0 and NPositive > 0) {
						// Positive weights
						Matrix[i][j] = PascalsDistribution(m, p, rand());
						NPositive--;
					}
					else {
						// Negative weights
						Matrix[i][j] = (-1) * PascalsDistribution(m, p, rand());
						NNegative--;
					}
					c++;
				}
			}
		}
	}
	// If weights can only be positive values
	else {
		// Go through every row
		for (int i = 0; i < Matrix.size(); i++) {
			// Go through every column
			for (int j = 0; j < Matrix[i].size(); j++) {
				// If vertices are adjacent
				if (AdjMatrix[i][j] == 1) {
					Matrix[i][j] = PascalsDistribution(m, p, rand());
				}
			}
		}
	}

}

string RouteAlgs();
string MSTAlgs();
void ProcessRouteAlg(OrientedGraph* OrGraph, string Alg);