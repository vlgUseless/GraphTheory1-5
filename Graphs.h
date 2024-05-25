#pragma once
#include <vector>
#include <map>
#include <fstream>
using namespace std;

class AbstractGraph
{
protected:
	unsigned int Nvertex; // Количество вершин
	unsigned int Nedge; // Количество рёбер

	vector<vector<int>> AdjMatrix; // Матрица смежности
	vector<vector<double>> WeightMatrix; // Матрица весов
	vector<vector<int>> ReachabilityMatrix; // Матрица достижимости
	vector<vector<int>> AdjList; // Список смежности
	vector<int> VertexDegreeList; // Список степеней вершин
	vector<vector<double>> CapacityMatrix; // Матрица пропускных способностей
	vector<vector<int>> KirchhoffMatrix; // Матрица Кирхгофа

	vector<int> Parent;
	vector<int> Rank;
	vector<int> Component;
	vector<pair<double, pair<int, int>>> EdgesList; // Список рёбер и их весов. Первый элемент вес, второй - пара вершин, образующих ребро
	vector<pair<double, pair<int, int>>> MST; // Остов графа
	vector<vector<int>> MSTmatrix;

	multimap<int, pair<int,double>> TreeView;
	
public:
	AbstractGraph(unsigned int _Nvertex = 0, unsigned int _Nedge = 0);
	AbstractGraph(const AbstractGraph& AbGraph);
	vector<vector<int>> getAdjMatrix() const;
	vector<vector<double>> getWeightMatrix() const;
	vector<vector<int>> getReachabilityMatrix() const;
	vector<vector<int>> getAdjList() const;
	vector<int> getVertexDegreeList() const;
	unsigned int getNvertex() const;
	unsigned int getNedge() const;
	vector<vector<int>> getMSTmatrix() const;

	void inputAdjMatrix();

	bool Dijkstra(vector<vector<double>> Matrix, int startV, int endV, vector<int>& route,  bool PrintRoute);
	void FloydWarshall(int startV, int endV);
	bool BFS(vector<vector<double>> rGraph, int s, int t, vector<int>& parent);
	bool MaxRoute(vector<vector<double>> Matrix, int startV, int endV, vector<int>& route, bool PrintRoute);
	void DFS(int startV, int endV, vector<bool>& visited, vector<int>& route);
	double SpanningTrees();

	void Kruskal();
	void Boruvka(bool print);
	int FindSet(int i, int& iter);
	void UnionSets(int x, int y, int& iter);

	vector<pair<int, double>> PruferEncode();
	multimap<int, pair<int, double>> PruferDecode(const vector<pair<int, double>>& Encoded);
};


class OrientedGraph : public AbstractGraph
{
public:
	OrientedGraph(unsigned int _Nvertex = 0, unsigned int _Nedge = 0);
	int FordFulkerson(int startV, int endV);
	double MinCostFlow();
};

class UnorientedGraph : public AbstractGraph
{
	vector<bool> visited;
	bool hasCycle;
public:
	UnorientedGraph(OrientedGraph* OrGraph);
	UnorientedGraph(const UnorientedGraph& graph);
	bool CheckEulerian();
	bool CheckHamiltonian();

	void MakeEulerian();
	vector<int> Fleury();

	void MakeHamiltonian();
	void HamiltonianCycles(ofstream& out);
	void FindHamCycle(vector<vector<int>> graph, int pos, vector<int> path, ofstream& out, double& ans);

};

vector<vector<double>> ShimbelsMethod(const vector<vector<double>>& WeightMatrix, unsigned int NumberOfEdges, const bool MaxRoute);
double ShimbelsRoute(const vector<vector<double>>& Matrix, const bool MaxRoute);

bool BFS(vector<vector<double>> rGraph, int s, int t, vector<int>& parent);

bool isSafe(int v, vector<vector<int>>& graph, vector<int>& path, int pos);