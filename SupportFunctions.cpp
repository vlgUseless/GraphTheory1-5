#include "SupportFunctions.h"
#include "Graphs.h"
#include <fstream>
using namespace std;

int InputVertices() {

	int vertices = 0;

	cout << "Enter the number of vertices of the graph: ";

	while (true) {
		cin >> vertices;

		if (cin.fail() or vertices < 2) {
			cin.clear();
			cin.ignore(numeric_limits<streamsize>::max(), '\n');
			cout << "Input must be a positive integer number from 2, try again: ";
		}

		else {
			break;
		}
	}
	return vertices;
}

int InputVertex(OrientedGraph* OrGraph, int numb) {
	int vertex = -1;

	cout << "Enter the vertex " << "#" << numb << ": ";

	while (true) {
		cin >> vertex;

		if (cin.fail() or vertex < 0 or vertex >= OrGraph->getAdjMatrix().size()) {
			cin.clear();
			cin.ignore(numeric_limits<streamsize>::max(), '\n');
			cout << "Input must be a positive integer number below " << OrGraph->getAdjMatrix().size() << ", try again: ";
		}

		else {
			break;
		}
	}

	return vertex;
}

int InputEdges(OrientedGraph* OrGraph) {
	int edges = 0;

	cout << "Enter the number of edges: ";

	while (true) {
		cin >> edges;

		if (cin.fail() or edges < 1 or edges > OrGraph->getAdjMatrix().size() - 1) {
			cin.clear();
			cin.ignore(numeric_limits<streamsize>::max(), '\n');
			cout << "Input must be a positive integer number below or equal " << OrGraph->getAdjMatrix().size() - 1 << ", try again: ";
		}

		else {
			break;
		}
	}

	return edges;
}

bool InputMinOrMax() {
	bool MinOrMax = 0;
	cout << "Enter 0 for the minimum route or 1 for the maximum route: ";
	while (true) {
		cin >> MinOrMax;

		if (cin.fail()) {
			cin.clear();
			cin.ignore(numeric_limits<streamsize>::max(), '\n');
			cout << "Input must be a boolean, try again: ";
		}

		else {
			break;
		}
	}
	return MinOrMax;
}


bool InputWeightPar() {
	bool WeightPar = 0;
	cout << "-Generating Weight Matrix-" << endl;
	cout << "[0] mixed weights (weights can be negative)" << endl;
	cout << "[1] positive weights" << endl << endl;
	while (true) {
		cin >> WeightPar;

		if (cin.fail()) {
			cin.clear();
			cin.ignore(numeric_limits<streamsize>::max(), '\n');
			cout << "Input must be a boolean, try again: ";
		}

		else {
			break;
		}
	}
	return WeightPar;
}
/*
*/
vector<vector<int>> MultiplyMatrix(const vector<vector<int>>& Matrix1, const vector<vector<int>>& Matrix2) {
	// If Matrixes are quadratic
	if (Matrix1.size() == Matrix2.size() and Matrix1[0].size() == Matrix2[0].size()) {
		int sz = Matrix1.size(); // Количество вершин
		vector<vector<int>> Result(sz, vector<int>(sz, 0));
		for (int i = 0; i < sz; i++) {
			for (int j = 0; j < sz; j++) {
				for (int k = 0; k < sz; k++) {
					Result[i][j] += (Matrix1[i][k] * Matrix2[k][j]);
				}

			}
		}
		return Result;
	}
}


void PrintList2d(const vector<vector<int>> List, string ListName) {
	cout << "---" << ListName << "---" << endl;
	for (int i = 0; i < List.size(); i++) {
		cout << i << ": ";
		for (int j = 0; j < List[i].size(); j++) {
			cout << List[i][j] << " ";
		}
		if (List[i].size() == 0) {
			cout << "-";
		}
		cout << endl;
	}
	cout << endl;
}




// Ручное заполнение матрицы весов
void InputWeightMatrix(vector<vector<double>>& Matrix, const vector<vector<int>>& AdjMatrix) {
	cout << "Enter the values of Weight Matrix for each rows. They must be based on the completion of the Adjacency Matrix: " << endl;

	// Go through every row
	for (int i = 0; i < Matrix.size(); i++) {
		// Go through every column
		for (int j = 0; j < Matrix[i].size(); j++) {
			// If vertices are adjacent
			if (AdjMatrix[i][j] == 1) {
				cin >> Matrix[i][j];
			}
		}
		cout << endl;
	}
}




void UnionMatrix(vector<vector<int>>& Matrix1, const vector<vector<int>>& Matrix2) {
	if (Matrix1.size() == Matrix2.size() and Matrix1[0].size() == Matrix2[0].size()) {
		int sz = Matrix1.size();
		for (int i = 0; i < sz; i++) {
			for (int j = 0; j < sz; j++) {
				Matrix1[i][j] = max(Matrix1[i][j], Matrix2[i][j]);
			}
		}
	}
}

vector<vector<int>> AddMatrix(const vector<vector<int>>& Matrix1, const vector<vector<int>>& Matrix2) {
	int sz = Matrix1.size();
	vector<vector<int>> Result(sz, vector<int>(sz, 0));
	for (int i = 0; i < sz; i++) {
		for (int j = 0; j < sz; j++) {
			Result[i][j] = Matrix1[i][j] + Matrix2[i][j];
		}
	}
	return Result;
}

string RouteAlgs() {
	string Alg = "0";
	cout << "-Choose an algorithm for finding route-" << endl;
	cout << "[0] Dijkstra (shortest route) /only for positive weights/" << endl;
	cout << "[1] Floyd-Warshall (shortest route)" << endl;
	cout << "[2] Max route algorithm" << endl << endl;
	while (true) {
		cin >> Alg;

		if (Alg.size() == 1 and (Alg == "0" or Alg == "1" or Alg == "2")) {
			break;
		}

		else {
			cin.clear();
			cin.ignore(numeric_limits<streamsize>::max(), '\n');
			cout << "Input must be a number 0, 1 or 2, try again: ";
		}
	}
	return Alg;
}


//TODO: обработка startV > endV внутри самих алгоритмов
void ProcessRouteAlg(OrientedGraph* OrGraph, string Alg) {
	vector<int> route;
	// If Dijkstra algorithm
	if (Alg == "0") {
		cout << "--Dijkstra algorithm--" << endl;
		int n = 0;
		int startV = InputVertex(OrGraph, n++);
		int endV = InputVertex(OrGraph, n++);
		
		if (startV <= endV) {
			OrGraph->Dijkstra(OrGraph->getWeightMatrix(), startV, endV, route, 1);
		}
		else {
			cout << "There's no route from vertex " << startV << " and vertex " << endV << "." << endl;
		}
	}
	// If Floyd-Warshall algorithm
	else if (Alg == "1") {
		cout << "--Floyd-Warshall algorithm--" << endl;
		int n = 0;
		int startV = InputVertex(OrGraph, n++);
		int endV = InputVertex(OrGraph, n++);
		if (startV <= endV) {
			OrGraph->FloydWarshall(startV, endV);
		}
		else {
			cout << "There's no route from vertex " << startV << " and vertex " << endV  << "." << endl;
		}
	}
	else if (Alg == "2") {
		cout << "--Longest route algorithm--" << endl;
		int n = 0;
		int startV = InputVertex(OrGraph, n++);
		int endV = InputVertex(OrGraph, n++);
		if (startV <= endV) {
			OrGraph->MaxRoute(OrGraph->getWeightMatrix() , startV, endV, route, 1);
			//int sm = 0;
			//for (int i = 0; i < route.size() - 1; i++) {
			//	sm += OrGraph->getWeightMatrix()[route[i]][route[i + 1]];
			//}
			//cout << "The longest route is " << sm << ": ";
			//for (int i = 0; i < route.size(); i++) {
			//	if (i != route.size() - 1) {
			//		cout << route[i] << "->";
			//	}
			//	else {
			//		cout << route[i] << endl << endl;
			//	}
			//}
		}
		else {
			cout << "There's no route from vertex " << startV << " and vertex " << endV << "." << endl;
		}
		
	}
}
string MSTAlgs() {
	string Alg = "0";
	cout << "-Choose an algorithm for finding MST-" << endl;
	cout << "[0] Kruskal's algorithm" << endl;
	cout << "[1] Boruvka's algorithm" << endl << endl;
	while (true) {
		cin >> Alg;

		if (Alg.size() == 1 and (Alg == "0" or Alg == "1")) {
			break;
		}

		else {
			cin.clear();
			cin.ignore(numeric_limits<streamsize>::max(), '\n');
			cout << "Input must be a number 0 or 1, try again: ";
		}
	}
	return Alg;
}

void ProcessMST(OrientedGraph* OrGraph, string Alg) {
	if (Alg == "0") {
		cout << "--Kruskal's algorithm--" << endl;
		OrGraph->Kruskal();
	}
	if (Alg == "1") {
		cout << "--Boruvka's algorithm--" << endl;
		OrGraph->Boruvka(true);
	}
}


void InputProcessing(string Action, OrientedGraph*& OrGraph) {
	UnorientedGraph* UnorGraph = new UnorientedGraph(OrGraph);
	// Regenerate OrGraph
	if (Action == "1") {
		
		delete OrGraph;
		delete UnorGraph;

		int vertices = InputVertices();
		OrGraph = new OrientedGraph(vertices);
		UnorGraph = new UnorientedGraph(OrGraph);
	}
	else if (Action == "2") {
		int edges = InputEdges(OrGraph);
		bool MinOrMax = InputMinOrMax();

		PrintMatrix(ShimbelsMethod(OrGraph->getWeightMatrix(), edges, MinOrMax), "Shimbels Matrix");

	}
	else if (Action == "3") {
		// Getting numbers of vertices [0; Nvertex-1]
		int numb = 1;
		int vertex1 = InputVertex(OrGraph, numb++);
		int vertex2 = InputVertex(OrGraph, numb++);

		PrintMatrix(OrGraph->getReachabilityMatrix(), "Reachability Matrix");

		// If route exists
		if (OrGraph->getReachabilityMatrix()[vertex1][vertex2] > 0) {
			cout << "Number of the routes: " << OrGraph->getReachabilityMatrix()[vertex1][vertex2] << endl << endl;
		}
		else {
			cout << "There's no route between vertices " << vertex1 << " and " << vertex2 << "." << endl;
		}
	}
	else if (Action == "4") {
		ProcessRouteAlg(OrGraph, RouteAlgs());
	}
	else if (Action == "5") {
		int numb = 1;
		int vertex1 = InputVertex(OrGraph, numb++);
		int vertex2 = InputVertex(OrGraph, numb++);
		cout << "Maximum flow is: " << OrGraph->FordFulkerson(vertex1, vertex2) << endl << endl;
	}
	else if (Action == "6") {
		cout << "---Minimal cost flow---" << endl;
		OrGraph->MinCostFlow();
		
	}
	else if (Action == "7") {
		cout << "The number of spanning trees of the graph: " << OrGraph->SpanningTrees() << endl << endl;
	}
	else if (Action == "8") {
		ProcessMST(OrGraph, MSTAlgs());
	}
	else if (Action == "9") {
		// Filling TreeView field
		OrGraph->Boruvka(false);

		vector<pair<int, double>> Encoded = OrGraph->PruferEncode();
		cout << "---Encoded by Prufer's code tree---" << endl;
		cout << "Code - Weight" << endl;
		for (int i = 0; i < Encoded.size(); i++) {
			cout << Encoded[i].first << " - " << Encoded[i].second << endl;
		}

		cout << "---Decoded by Prufer's code tree---" << endl;
		multimap<int, pair<int, double>> Decoded = OrGraph->PruferDecode(Encoded);
		for (const auto& el : Decoded) {
			cout << "(" << el.first << ", " << el.second.first << ") - " << el.second.second << endl;
		}

		//vector<vector<int>> MSTmatrix = OrGraph->getMSTmatrix();
		//for (int i = 0; i < MSTmatrix.size(); i++) {
		//	for (int j = 0; j < MSTmatrix.size(); j++) {
		//		if (MSTmatrix[i][j])
		//	}

		//}
	}
	else if (Action == "10") {
		cout << "---Checking unoriented graph---" << endl;
		cout << "Eulerian graph: " << boolalpha << UnorGraph->CheckEulerian() << endl;
		cout << "Hamiltonian graph: " << boolalpha << UnorGraph->CheckHamiltonian() << endl << endl;
	}
	// Eulerian graph
	else if (Action == "11") {
		UnorientedGraph* FictGraph = new UnorientedGraph(*UnorGraph);

		if (FictGraph->getNvertex() == 2) {
			cout << "Graph can't be Eulerian because of 2 vertices." << endl << endl;
			return;
		}

		cout << "-----Original graph-----" << endl << endl;
		PrintMatrix(UnorGraph->getAdjMatrix(), "Adjacency Matrix");
		PrintMatrix(UnorGraph->getWeightMatrix(), "Weight matrix");
		PrintList(UnorGraph->getVertexDegreeList(), "List of vertex degrees");
		PrintList2d(UnorGraph->getAdjList(), "Adjacency List");

		if (UnorGraph->CheckEulerian()) {
			cout << "The graph is already Eulerian" << endl << endl;
		}
		else {

			FictGraph->MakeEulerian();

			cout << "-----New Eulerian graph-----" << endl << endl;
			PrintMatrix(FictGraph->getAdjMatrix(), "Adjacency Matrix");
			PrintMatrix(FictGraph->getWeightMatrix(), "Weight matrix");
			PrintList(FictGraph->getVertexDegreeList(), "List of vertex degrees");
			PrintList2d(FictGraph->getAdjList(), "Adjacency List");
			cout << "Is Eulerian: " << boolalpha << FictGraph->CheckEulerian() << endl << endl;
		}

		vector<int> route = FictGraph->Fleury();
		double w = 0;
		cout << "Eulerian path: <";
		for (int i = 0; i < route.size(); i++) {
			if (i != route.size() - 1) {
				w += FictGraph->getWeightMatrix()[route[i]][route[i + 1]];
				cout << route[i] << ";";
			}
			else {
				cout << route[i] << ">, weight = " << w << endl << endl;
			}
		}
	}
	// Hamiltonian graph
	else if (Action == "12") {
		UnorientedGraph* FictGraph = new UnorientedGraph(*UnorGraph);

		if (FictGraph->getNvertex() == 2) {
			cout << "Graph can't be Hamiltonian because of 2 vertices." << endl << endl;
			return;
		}

		cout << "-----Original graph-----" << endl << endl;
		PrintMatrix(UnorGraph->getAdjMatrix(), "Adjacency Matrix");
		PrintMatrix(UnorGraph->getWeightMatrix(), "Weight Matrix");
		PrintList(UnorGraph->getVertexDegreeList(), "List of vertex degrees");
		PrintList2d(UnorGraph->getAdjList(), "Adjacency List");

		if (UnorGraph->CheckHamiltonian()) {
			cout << "The graph is already Hamiltonian" << endl << endl;
		}
		else {
			FictGraph->MakeHamiltonian();

			cout << "-----New Hamiltonian graph-----" << endl << endl;
			PrintMatrix(FictGraph->getAdjMatrix(), "Adjacency Matrix");
			PrintMatrix(FictGraph->getWeightMatrix(), "Weight Matrix");
			PrintList(FictGraph->getVertexDegreeList(), "List of vertex degrees");
			PrintList2d(FictGraph->getAdjList(), "Adjacency List");
			cout << "Is Hamiltonian: " << boolalpha << FictGraph->CheckHamiltonian() << endl << endl;

		}
		ofstream out("Hamiltonian cycles.txt");
		//FictGraph->inputAdjMatrix();
		FictGraph->HamiltonianCycles(out);
	}
}

void Menu() {
	int vertices = InputVertices();
	OrientedGraph* OrGraph = new OrientedGraph(vertices);
	while (true)
	{
		cout << endl;
		cout << "Select an action:" << endl;
		cout << "(1) Regenerate an oriented connected acyclic graph" << endl;
		cout << "(2) Use Shimbel's method on generated graph" << endl;
		cout << "(3) Check the feasibility of the route" << endl;
		cout << "(4) Find the route" << endl;
		cout << "(5) Find the maximum flow by Ford-Fulkerson's algorithm" << endl;
		cout << "(6) Find the minimal cost for the 2/3 of the max flow" << endl;
		cout << "(7) Find the number of spanning trees of the graph" << endl;
		cout << "(8) Find the minimal weight spanning tree of the graph" << endl;
		cout << "(9) Encode and decode MST with Prufer's code" << endl;
		cout << "(10) Check if generated unoriented graph is Euler or Hamiltonian" << endl;
		cout << "(11) Modify the graph to be Eulerian and find Eulerian cycles" << endl;
		cout << "(12) Modify the graph to be Hamiltonian and find Hamiltonian cycles" << endl;
		cout << "(0) Exit " << endl << endl;

		string InputAction;
		cin >> InputAction;
		cin.clear();
		cin.ignore(numeric_limits<streamsize>::max(), '\n');

		if (InputAction == "0")
		{
			break;
		}
		cout << endl;

		InputProcessing(InputAction, OrGraph);

		cout << endl;
	}
}

// Проверка метода Шимбелла на любой матрице весов
void TestShimbelsMethod() {
	vector<vector<double>> Weight = InputMatrix<double>();

	int edges = 0;
	cout << "Input the number of edges: ";
	cin >> edges;
	bool MinOrMax = InputMinOrMax();

	PrintMatrix(ShimbelsMethod(Weight, edges, MinOrMax), "Shimbel's Matrix");
}

double GaussDet(vector<vector<double>> m)
{
	int n = m.size();
	double det = 1;
	for (int i = 0; i < n; ++i)
	{
		double mx = fabs(m[i][i]);
		int idx = i;
		for (int j = i + 1; j < n; ++j)
			if (mx < fabs(m[i][j])) mx = fabs(m[i][idx = j]);
		if (idx != i)
		{
			for (int j = i; j < n; ++j)
			{
				double t = m[j][i];
				m[j][i] = m[j][idx];
				m[j][idx] = t;
			}
			det = -det;
		}
		for (int k = i + 1; k < n; ++k)
		{
			double t = m[k][i] / m[i][i];

			for (int j = i; j < n; ++j)
				m[k][j] -= m[i][j] * t;
		}
	}
	for (int i = 0; i < n; ++i) det *= m[i][i];
	return det;
}