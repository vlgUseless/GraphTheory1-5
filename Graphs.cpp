#include "Graphs.h"
#include "SupportFunctions.h"
#include <random>
#include <tuple>
#include <queue>
#include <iostream>
#include <functional>
#include <stack>
using namespace std;


AbstractGraph::AbstractGraph(unsigned int _Nvertex, unsigned int _Nedge): Nvertex(_Nvertex), Nedge(_Nedge)
{
	AdjMatrix.resize(Nvertex, vector<int>(Nvertex, 0));
	WeightMatrix.resize(Nvertex, vector<double>(Nvertex, INFINITY));
	ReachabilityMatrix.resize(Nvertex, vector<int>(Nvertex, 0));
	AdjList.resize(Nvertex, vector<int>());
	AdjMatrix.resize(Nvertex, vector<int>(Nvertex, 0));
	VertexDegreeList.resize(Nvertex);
	KirchhoffMatrix.resize(Nvertex, vector<int>(Nvertex, 0));
	CapacityMatrix.resize(Nvertex, vector<double>(Nvertex, 0));
	Parent.resize(Nvertex, -1);
	Rank.resize(Nvertex, -1);
	EdgesList.resize(Nedge);
	Component.resize(Nvertex, -1);
	MSTmatrix.resize(Nvertex, vector<int>(Nvertex, 0));
}

AbstractGraph::AbstractGraph(const AbstractGraph & AbGraph) {
	Nvertex = AbGraph.Nvertex;
	Nedge = AbGraph.Nedge;
	AdjMatrix = AbGraph.AdjMatrix;
	WeightMatrix = AbGraph.WeightMatrix;
	ReachabilityMatrix = AbGraph.ReachabilityMatrix;
	AdjList = AbGraph.AdjList;
	VertexDegreeList = AbGraph.VertexDegreeList;
	CapacityMatrix = AbGraph.CapacityMatrix;
	KirchhoffMatrix = AbGraph.KirchhoffMatrix;
	EdgesList = AbGraph.EdgesList;
	MST = AbGraph.MST;
	TreeView = AbGraph.TreeView;
}

vector<vector<int>> AbstractGraph::getAdjMatrix() const {return AdjMatrix;}
vector<vector<double>> AbstractGraph::getWeightMatrix() const { return WeightMatrix; }
vector<vector<int>> AbstractGraph::getReachabilityMatrix() const { return ReachabilityMatrix; }
vector<vector<int>> AbstractGraph::getAdjList() const { return AdjList; }
unsigned int AbstractGraph::getNvertex() const { return Nvertex; }
unsigned int AbstractGraph::getNedge() const { return Nedge; }
vector<int> AbstractGraph::getVertexDegreeList() const {return VertexDegreeList;}
vector<vector<int>> AbstractGraph::getMSTmatrix() const { return MSTmatrix; }

void AbstractGraph::inputAdjMatrix() {
	AdjMatrix = InputMatrix<int>();
}


// Генерация случайного графа
OrientedGraph::OrientedGraph(unsigned int _Nvertex, unsigned int _Nedge) : AbstractGraph(_Nvertex, _Nedge)
{
	// Параметры генерации случайных чисел
	srand(time(0));
	int m = Nvertex/2;
	double p = 0.9;
	
	// Генерация количества рёбер, исходящих из каждой вершины / Заполнение вектора Parent
	int sm = 0;
	int i = 0;
	vector<int> tmpVDL(Nvertex, 0);
	while (i < Nvertex) {
		int gnumber = PascalsDistribution(m, p, rand());
		if (gnumber <= Nvertex - 1) {
			tmpVDL[i] = gnumber;
			Parent[i] = i;
			Rank[i] = 0;
			i++;
		}
	}

	
	PrintList<int>(tmpVDL, "Generated Vertices Degrees List");

	// Заполнение матрицы и списка смежности
	for (int i = 0; i < Nvertex; i++) {
		for (int j = i + 1; j < Nvertex; j++) {
			if (i == j) { AdjMatrix[i][j] = 0; continue; }
			if ((tmpVDL[i]--) > 0) {
				Nedge++;
				AdjMatrix[i][j] = 1;
				AdjList[i].push_back(j);
			}
		}
	}
	for (int i = 0; i < Nvertex; i++) {
		for (int j = 0; j < Nvertex; j++) {
			VertexDegreeList[i] += AdjMatrix[i][j] + AdjMatrix[j][i];
		}
	}

	PrintMatrix(getAdjMatrix(), "Adjacency Matrix");
	PrintList2d(getAdjList(), "Adjacency List");


	bool WeightPar = InputWeightPar();
	GenerateMatrix(WeightMatrix, AdjMatrix, Nedge, 7, 0.3, WeightPar);
	for (int i = 0; i < Nvertex; i++) {
		for (int j = 0; j < Nvertex; j++) {
			if (AdjMatrix[i][j] > 0) {
				EdgesList.push_back(make_pair(WeightMatrix[i][j], make_pair(i, j)));
			}
		}
	}
	PrintMatrix(WeightMatrix, "Weight Matrix");
	
	// Заполнение матрицы достижимости
	vector<vector<int>> TempMatrix = AdjMatrix;
	for (int i = 0; i < Nvertex; i++) {
		// Заполняем главную диагональ
		ReachabilityMatrix[i][i] = 1;
	}
	for (int i = 0; i < Nvertex - 1; i++) {
		ReachabilityMatrix = AddMatrix(ReachabilityMatrix, TempMatrix);
		TempMatrix = MultiplyMatrix(TempMatrix, AdjMatrix);
	}

	GenerateMatrix(CapacityMatrix, AdjMatrix, Nedge, 5, 0.5, 1);
	PrintMatrix(CapacityMatrix, "Capacity Matrix");

	// Filling Kirchhoff's Matrix
	for (int i = 0; i < Nvertex; i++) {
		for (int j = 0; j < Nvertex; j++) {
			if (i == j) {
				KirchhoffMatrix[i][j] = VertexDegreeList[i];
			}
			else {
				KirchhoffMatrix[i][j] = -AdjMatrix[i][j];
			}
		}
	}
	PrintMatrix(KirchhoffMatrix, "Kirchhoff's Matrix");
}

UnorientedGraph::UnorientedGraph(OrientedGraph* OrGraph): AbstractGraph(OrGraph->getNvertex(), OrGraph->getNedge()) {
	hasCycle = false;
	visited.resize(Nvertex, 0);
	VertexDegreeList = OrGraph->getVertexDegreeList();
	AdjMatrix = OrGraph->getAdjMatrix();
	WeightMatrix = OrGraph->getWeightMatrix();

	// Дополняем граф до неориентированного
	for (int i = 0; i < Nvertex; i++) {
		for (int j = i + 1; j < Nvertex; j++) {
			AdjMatrix[j][i] = AdjMatrix[i][j];
			WeightMatrix[j][i] = WeightMatrix[i][j];
		}
	}
	for (int i = 0; i < Nvertex; i++) {
		for (int j = 0; j < Nvertex; j++) {
			if (AdjMatrix[i][j]) {
				AdjList[i].push_back(j);
			}
			
		}
	}
}

UnorientedGraph::UnorientedGraph(const UnorientedGraph& graph){
	hasCycle = false;
	visited.resize(Nvertex, 0);
	*this = graph;
}

// Args:
// - MinOrMax:
// -- 0 if Min, 1 if Max
vector<vector<double>> ShimbelsMethod(const vector<vector<double>>& WeightMatrix, unsigned int NumberOfEdges, const bool MaxRoute) {
	vector<vector<double>> ShimbelsMatrix = WeightMatrix;
	int sz = WeightMatrix.size();
	
	// Doing the multiplying of matrices (NumberOfEdges - 1) times
	for (int edges = 0; edges < NumberOfEdges - 1; edges++) {
		vector<vector<double>> SupportMatrix = ShimbelsMatrix;
		// Going through rows of ShimbelsMatrix
		for (int i = 0; i < sz; i++) {
			// Going through columns of WeightMatrix
			for (int j = 0; j < sz; j++) {
				
				vector<double> Elements;
				// Going through elements of rows and columns of Shimbels(Max/Min)Matrix and WeightMatrix
				for (int k = 0; k < sz; k++) {
					if (SupportMatrix[i][k] != INFINITY and WeightMatrix[k][j] != INFINITY) {
						Elements.push_back(SupportMatrix[i][k] + WeightMatrix[k][j]);
					}
				}
				if (Elements.empty()) {
					Elements.push_back(INFINITY);
				}

				if (MaxRoute == true) {
					ShimbelsMatrix[i][j] = *max_element(Elements.begin(), Elements.end());
				}
				else {
					ShimbelsMatrix[i][j] = *min_element(Elements.begin(), Elements.end());
				}
				
			}
		}
	}
	return ShimbelsMatrix;
}

double ShimbelsRoute(const vector<vector<double>>& Matrix, const bool MaxRoute) {
	double mn = 2 << 15;
	double mx = -(2 << 15);
	for (int i = 0; i < Matrix.size(); i++) {
		for (int j = 0; j < Matrix.size(); j++) {
			if (Matrix[i][j] != INFINITY) {
				mn = min(mn,Matrix[i][j]);
				mx = max(mx,Matrix[i][j]);
			}
		}
	}
	if (MaxRoute == true) {
		return mx;
	}
	return mn;
}



bool AbstractGraph::Dijkstra(vector<vector<double>> Matrix, int startV, int endV, vector<int>& route, bool PrintRoute) {
	int N = Matrix.size();
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			if (Matrix[i][j] < 0) {
				cout << "This implementation of the Dijkstra's algorithm doesn't work with the negative weights." << endl << endl;
				return 0;
			}
		}
	}

	vector<double> T(N, INFINITY); // Вектор со значениями кратчайших путей
	vector<int> X(N, 0); // Все вершины изначально не отмечены
	vector<int> H(N, -1); // Вектор вершин, предшествующих i на кратчайшем пути
	T[startV] = 0; // Кратчайший путь изначально имеет длину 0
	X[startV] = 1; // И он известен
	int curV = startV; // Текущая вершина

	int iter = 0;
	// Обновление меток
	while (true) {
		// Перебираем все вершины
		for (int u = 0; u < N; u++) {
			// Если в данной вершине мы ещё не были И в неё вообще можно попасть из текущей точки И
			// путь из текущей точки в неё короче прошлого пути
			if (X[u] == 0 and Matrix[curV][u] != INFINITY and (T[curV] + Matrix[curV][u]) < T[u]) {
				T[u] = T[curV] + Matrix[curV][u];
				H[u] = curV;
			}
		}
		double m = INFINITY;
		curV = -1;

		// Ищем кратчайший путь на текущей итерации
		for (int u = 0; u < N; u++) {
			// Если в данной вершине мы ещё не были И путь в неё меньше предыдущего
			if (X[u] == 0 and T[u] < m) {
				// Переходим в данную вершину
				curV = u;
				// Обновляем кратчайший путь
				m = T[u];
			}
		}

		if (curV == -1) {
			if (PrintRoute) {
				cout << "There's no route between vertices " << startV << " and " << endV << "." << endl;
			}
			
			return 0;
		}
		if (curV == endV) {
			route.push_back(endV);
			// Пока не дойдём до начальной вершины
			while (route[route.size()-1] != startV) {
				route.push_back(H[route[route.size() - 1]]);
			}
			break;

		}
		X[curV] = 1;
		iter++;
	}
	sort(route.begin(), route.end());
	if (PrintRoute) {
		cout << "Shortest route is " << T[curV] << ": ";
		for (int i = 0; i < route.size(); i++) {
			if (i != route.size() - 1) {
				cout << route[i] << "->";
			}
			else {
				cout << route[i] << endl;
			}
		}
		cout << "Number of iterations: " << iter * 2 * Nvertex << endl << endl;

		PrintList<double>(T, "Vector of shortest routes");
	}
	return 1;
}

void AbstractGraph::FloydWarshall(int startV, int endV) {
	vector<vector<double>> T = WeightMatrix; // Матрица длин путей
	vector<vector<int>> H = AdjMatrix; // Матрица самих путей

	int sz = AdjMatrix.size();
	for (int i = 0; i < sz; i++) {
		for (int j = 0; j < sz; j++) {
			if (H[i][j] > 0) {
				H[i][j] = j; // Есть дуга из i в j
			}
		}
	}

	for (int i = 0; i < sz; i++) {
		for (int j = 0; j < sz; j++) {
			for (int k = 0; k < sz; k++) {
				if (i != j and T[j][i] != INFINITY and i != k and T[i][k] != INFINITY and (T[j][k] == INFINITY or T[j][k] > T[j][i] + T[i][k])) {
					H[j][k] = H[j][i]; // Запомнить новый путь
					T[j][k] = T[j][i] + T[i][k]; // и его длину
				}
			}
			
		}
		
		for (int j = 0; j < sz; j++) {
			if (T[i][j] < 0) {
				//cout << "There is no solution." << endl;
				break; // Нет решения
			}
		}
	}

	// Извлекаем путь
	int w = startV;
	cout << "Shortest route is " << T[startV][endV] << ": ";
	while (w != endV) {
		cout << w << "->";
		w = H[w][endV];
	}
	cout << endV << endl;
	cout << "Number of iterations: " << sz*sz + sz*2*sz*sz<< endl << endl;
	PrintMatrix(T, "Distance matrix");
}


void AbstractGraph::DFS(int startV, int endV, vector<bool>& visited, vector<int>& route) {
	visited[startV] = true;
	route.push_back(startV);

	for (int i = startV; i != endV + 1; i++)
		if (!visited[i]) {
			DFS(i, endV, visited, route);
		}
}

bool AbstractGraph::MaxRoute(vector<vector<double>> Matrix, int startV, int endV, vector<int>& route, bool PrintRoute) {
	int N = Matrix.size();
	vector<double> dist(N, -INFINITY); // Вектор со значениями максимальных путей
	vector<int> parent(N, -1); // Вектор предшествующих вершин
	dist[startV] = 0; // Максимальный путь до начальной вершины - 0

	// Топологическая сортировка вершин
	vector<int> topOrder;
	vector<bool> visited(N, false);
	function<void(int)> topologicalSort = [&](int v) {
		visited[v] = true;
		for (int i = 0; i < N; ++i) {
			if (Matrix[v][i] != INFINITY && !visited[i]) {
				topologicalSort(i);
			}
		}
		topOrder.push_back(v);
	};

	for (int i = 0; i < N; ++i) {
		if (!visited[i]) {
			topologicalSort(i);
		}
	}
	reverse(topOrder.begin(), topOrder.end());


	for (int u : topOrder) {
		if (dist[u] != -INFINITY) {
			for (int v = 0; v < N; ++v) {
				if (Matrix[u][v] != INFINITY) {
					if (dist[u] + Matrix[u][v] > dist[v]) {
						dist[v] = dist[u] + Matrix[u][v];
						parent[v] = u;
					}
				}
			}
		}
	}

	// Формирование пути
	if (dist[endV] == -INFINITY) {
		if (PrintRoute) {
			cout << "There's no route between vertices " << startV << " and " << endV << "." << endl;
		}
		return false;
	}
	else {
		for (int v = endV; v != -1; v = parent[v]) {
			route.push_back(v);
		}
		reverse(route.begin(), route.end());

		if (PrintRoute) {
			cout << "Longest route is " << dist[endV] << ": ";
			for (size_t i = 0; i < route.size(); ++i) {
				if (i != route.size() - 1) {
					cout << route[i] << "->";
				}
				else {
					cout << route[i] << endl;
				}
			}
		}
		return true;
	}
}


double AbstractGraph::SpanningTrees()
{
	// Вычисляем матрицу с алгебраическим дополнением
	vector<vector<double>> M(Nvertex - 1,vector<double>(Nvertex-1,0));
	for (int i = 1; i < Nvertex; i++) {
		for (int j = 1; j < Nvertex; j++) {
			M[i - 1][j - 1] = KirchhoffMatrix[i][j];
		}
	}
	return GaussDet(M);
}

int AbstractGraph::FindSet(int i, int& iter) {
	iter++;
	if (i == Parent[i]) {
		return i;
	}
	return FindSet(Parent[i], iter);
}

void AbstractGraph::UnionSets(int x, int y, int& iter) {
	int rx = FindSet(x, iter);
	int ry = FindSet(y, iter);
	if (rx == ry) { return; }
	if (Rank[rx] > Rank[ry]) {
		Parent[ry] = rx;
	}
	else {
		Parent[rx] = ry;
		if (Rank[rx] == Rank[ry]) {
			Rank[ry] = Rank[rx] + 1;
		}
	}
}


void AbstractGraph::Kruskal()
{	
	int iter = 0;

	MST.clear();
	sort(EdgesList.begin(), EdgesList.end());
	iter += EdgesList.size() * log(EdgesList.size());

	cout << endl << "---Original graph---" << endl;
	cout << "Weight - Edge" << endl;
	for (const auto el : EdgesList) {
		cout << setw(6) << el.first << " - (" << el.second.first << "," << el.second.second << ")" << endl;
	}
	cout << endl;

	vector<int> Parent;
	vector<int> Rank;
	double weight = 0;

	int k = 0;
	for (int i = 0; i < EdgesList.size(); i++) {
		if (FindSet(EdgesList[i].second.first, iter) != FindSet(EdgesList[i].second.second, iter)) {
			weight += EdgesList[i].first;
			MST.push_back(EdgesList[i]);
			UnionSets(EdgesList[i].second.first, EdgesList[i].second.second, iter);
		}
	}

	cout << "Minimal weight spanning tree is " << weight << ": ";
	for (int i = 0; i < MST.size(); i++) {
		if (i != MST.size() - 1) {
			cout << "(" << MST[i].second.first << ", " << MST[i].second.second << "), ";
		}
		else {
			cout << "(" << MST[i].second.first << ", " << MST[i].second.second << ")" << endl;
		}
	}
	cout << "Number of iterations: " << iter << endl << endl;
}

void AbstractGraph::Boruvka(bool print)
{
	int iter = 0;

	TreeView.clear();
	MST.clear();
	MST.resize(Nvertex, make_pair(-1, make_pair(-1, -1)));
	//for (int i = 0; i < MST.size(); i++) {
	//	cout << i << ": " << MST[i].first << " " << MST[i].second.first << " " << MST[i].second.second << endl;
	//}

	// Initially there are V different trees.
	// Finally there will be one tree that will be MST
	int numTrees = Nvertex;
	int MSTweight = 0;
	// Create V subsets with single elements
	for (int node = 0; node < Nvertex; node++) {
		Parent[node] = node;
		Rank[node] = 0;
		iter++;
	}

	// Keep combining components (or sets) until all
	// components are not combined into single MST
	while (numTrees > 1) {
		// Traverse through all edges and update
		// cheapest of every component
		for (int i = 0; i < EdgesList.size(); i++) {
			// Find components (or sets) of two corners
			// of current edge
			int u = EdgesList[i].second.first, v = EdgesList[i].second.second,
				w = EdgesList[i].first;
			int set1 = FindSet(u, iter),
				set2 = FindSet(v, iter);

			// If two corners of current edge belong to
			// same set, ignore current edge. Else check
			// if current edge is closer to previous
			// cheapest edges of set1 and set2
			if (set1 != set2) {
				if (MST[set1].first == -1 || MST[set1].first > w) {
					MST[set1] = make_pair(w, make_pair(u, v));
				}
				if (MST[set2].first == -1 || MST[set2].first > w) {
					MST[set2] = make_pair(w, make_pair(u, v));
				}
			}
		}

		// Consider the above picked cheapest edges and
		// add them to MST
		for (int node = 0; node < Nvertex; node++) {
			// Check if cheapest for current set exists
			if (MST[node].first != -1) {
				int u = MST[node].second.first,
					v = MST[node].second.second,
					w = MST[node].first;
				int set1 = FindSet(u, iter),
					set2 = FindSet(v, iter);
				if (set1 != set2) {
					MSTweight += w;
					UnionSets(set1, set2, iter);
					if (print) {
						cout << "Edge " << u << "-" << v << " with weight " << w << " included in MST\n\n";
					}
					
					MSTmatrix[u][v] = WeightMatrix[u][v];
					TreeView.insert(make_pair(u, make_pair(v, w)));
					numTrees--;
				}
			}
		}

		for (int node = 0; node < Nvertex; node++) {
			// reset cheapest array
			MST[node].first = -1;
			iter++;
		}
	}
	if (print) {
		PrintMatrix(MSTmatrix, "MST Weight Matrix");
		cout << "Weight of MST is " << MSTweight << endl;
		cout << "Number of iterations: " << iter << endl << endl;
	}
}

bool AbstractGraph::BFS(vector<vector<double>> rGraph, int s, int t, vector<int>& parent)
{
	// Create a visited array and mark all vertices as not
	// visited
	int V = rGraph.size();
	vector<bool> visited(V, 0);

	// Create a queue, enqueue source vertex and mark source
	// vertex as visited
	queue<int> q;
	q.push(s);
	visited[s] = true;
	parent[s] = -1;

	// Standard BFS Loop
	while (!q.empty()) {
		int u = q.front();
		q.pop();

		for (int v = 0; v < V; v++) {
			if (visited[v] == false && rGraph[u][v] > 0) {
				// If we find a connection to the sink node,
				// then there is no point in BFS anymore We
				// just have to set its parent and can return
				// true
				if (v == t) {
					parent[v] = u;
					return true;
				}
				q.push(v);
				parent[v] = u;
				visited[v] = true;
			}
		}
	}

	// We didn't reach sink in BFS starting from source, so
	// return false
	return false;
}

int OrientedGraph::FordFulkerson(int startV, int endV)
{
	int u, v;
	int MaxFlow = 0;
	vector<vector<double>> rGraph = CapacityMatrix;

	vector<int> parent(Nvertex, 0); 

	// Augment the flow while there is path from source to
	// sink
	while (BFS(rGraph, startV, endV, parent)) {
		// Find minimum residual capacity of the edges along
		// the path filled by BFS. Or we can say find the
		// maximum flow through the path found.
		double path_flow = INFINITY;
		for (v = endV; v != startV; v = parent[v]) {
			u = parent[v];
			path_flow = min(path_flow, rGraph[u][v]);
		}

		// update residual capacities of the edges and
		// reverse edges along the path
		for (v = endV; v != startV; v = parent[v]) {
			u = parent[v];
			rGraph[u][v] -= path_flow;
			rGraph[v][u] += path_flow;
		}

		// Add path flow to overall flow
		MaxFlow += path_flow;
	}

	// Return the overall flow
	return MaxFlow;
}

double OrientedGraph::MinCostFlow()
{
	int mxflow = FordFulkerson(0, Nvertex - 1);
	int flow = (double)2 / 3 * mxflow; //TODO: передавать вершины и сразу считать фалкерсона?
	cout << "Max flow: " << mxflow << endl;
	cout << "Desired flow: " << flow << endl << endl;
	if (flow == 0) {
		return 0;
	}
	double CurFlow = 0;
	vector<vector<int>> all_path;
	vector<int> all_fi;

	vector<vector<double>> tmpCostMatrix = WeightMatrix;
	vector<vector<double>> tmpCapMatrix = CapacityMatrix; // Остаточная сеть/матрица

	vector<int> route;

	// TODO: перевернуть вектор пути и снова поправить индексы
	while (Dijkstra(tmpCostMatrix, 0, tmpCostMatrix.size() - 1, route, 0))
	{
		all_path.push_back(route);
		vector<int> tmp;
		for (int i = 0; i < route.size() - 1; i++)
		{
			tmp.push_back(tmpCapMatrix[route[i]][route[i + 1]]); // Добавляем пропускные способности каждого ребра пути
		}

		// Прибавляем к текущему потоку минимальную пропускную способность пути
		double mn = *min_element(begin(tmp), end(tmp));

		// Если перепрыгнули нужный поток, то добавляем в вектор разность и прекращаем поиск
		if (CurFlow + mn > flow)
		{
			all_fi.push_back(flow - CurFlow);
			break;
		}
		// Добавляем в вектор взятую пропускную способность
		CurFlow += mn;
		all_fi.push_back(mn);

		// Заполняем обратные рёбра остаточной сети как минимальный взятый поток
		for (int i = 0; i < route.size() - 1; i++)
		{
			tmpCapMatrix[route[i + 1]][route[i]] = mn;
		}


		for (int i = 0; i < route.size() - 1; i++)
		{
			// Для всех ненасыщенных дуг
			if (tmpCapMatrix[route[i]][route[i + 1]] != mn)
			{
				// Обновляем стоимость и поток
				tmpCostMatrix[route[i]][route[i + 1]] = WeightMatrix[route[i]][route[i + 1]];
				tmpCapMatrix[route[i]][route[i + 1]] = CapacityMatrix[route[i]][route[i + 1]] - mn;
			}
			// Для насыщенных дуг
			else
			{
				tmpCostMatrix[route[i]][route[i + 1]] = INFINITY;
				tmpCapMatrix[route[i]][route[i + 1]] = 0;
			}
		}
		route.clear();
	}

	cout << "-Counting minimal cost-" << endl;
	cout << setw(10) << "Edge" << setw(15) << "Current cost" << setw(20) << "Cost of the edge" << setw(15) << "Let in flow" << setw(15)<< "Current flow" << endl;
	int sumi = 0;
	int fl = 0;
	for (int i = 0; i < all_fi.size(); i++)
	{
		fl += all_fi[i];
		for (int j = 0; j < all_path[i].size() - 1; j++)
		{
			int w = WeightMatrix[all_path[i][j]][all_path[i][j + 1]] * all_fi[i];
			sumi += w;
			// Выводим данные
			//cout << "Current let in flow: " << all_fi[i] << endl;
			//cout << "General flow: " << sumi << endl;
			//PrintList(all_path[i], "Current route");
			if (all_fi[i] != 0) {
				cout << "(" << all_path[i][j] << ", " << all_path[i][j + 1] << ")" << setw(5) << " " << setw(15) << sumi << setw(20) << w \
					<< setw(15) << all_fi[i] << setw(15) << fl << endl;
			}
			
		}
	}
	cout << endl << "Total cost: " << sumi << endl;
	return sumi;
}

vector<pair<int, double>> AbstractGraph::PruferEncode() {
	vector<pair<int, double>> Encode;
	vector<int> DegreeList(Nvertex, 0);
	multimap<int, pair<int, double>> tempTree = TreeView;

	vector<int> deleted(Nvertex, 1);
	for (const auto& el : tempTree) {
		DegreeList[el.first]++;
		DegreeList[el.second.first]++;
	}

	cout << "--Current tree--" << endl;

	for (int i = 0; i < Nvertex-1; i++) {
		cout << "Iteration " << i << ":" << endl;
		for (const auto& el : tempTree) {
			cout << "(" << el.first << ", " << el.second.first << ") - " << el.second.second << endl;
		}


		// Ищем минимальную висячую вершину
		int v = Nvertex;
		cout << endl;
		for (int j = 0; j < Nvertex; j++) {
			//cout << j << ": " << DegreeList[j] << endl;
			if (DegreeList[j] == 1 and j < v) {
				v = j;
			}
		}
		cout << endl;

		// Кодируем вершину v
		for (auto it = tempTree.begin(); it != tempTree.end(); ++it) {
			if (it->second.first == v) {
				Encode.push_back(make_pair(it->first, it->second.second));
				DegreeList[v]--;
				DegreeList[it->first]--;
				tempTree.erase(it);
				break;
			}
			else if (it->first == v) {
				Encode.push_back(make_pair(it->second.first, it->second.second));
				DegreeList[v]--;
				DegreeList[it->second.first]--;
				tempTree.erase(it);
				break;
			}
		}
	}
	return Encode;
}



multimap<int, pair<int, double>> AbstractGraph::PruferDecode(const vector<pair<int, double>>& Encoded) {
	// Исходный граф в виде дерева
	multimap<int, pair<int, double>> Tree;
	vector<pair<int, double>> tmpPrufer = Encoded;
	// Список неиспользованных вершин
	vector<bool> B(Nvertex, 1);

	// Перебираем каждый из кодов
	for (int i = 0; i < Encoded.size(); i++) {
		// Ищем минимальную неиспользованную вершину, которая не встречается в коде Прюфера
		for (int j = 0; j < Nvertex; j++) {
			// Если вершина не использована
			if (B[j]) {
				bool flag = 1; // По умолчанию вершина подходит
				// Ищем вершину в коде Прюфера
				for (int k = i; k < tmpPrufer.size(); k++) {
					// Если текущая вершина есть в коде Прюфера, переходим к следующей
					if (Encoded[k].first == j) {
						flag = 0;
						break;
					}
				}

				// Если вершина подходит, она и есть минимальная
				if (flag) {
					Tree.insert(make_pair(Encoded[i].first, make_pair(j, Encoded[i].second)));
					B[j] = 0;
					break;
				}
			}

		}
	}

	return Tree;
}

//-----5 лабораторная-----
bool UnorientedGraph::CheckEulerian()
{
	if (Nvertex == 2) {
		return false;
	}
	bool check = true;
	for (int i = 0; i < VertexDegreeList.size(); i++) {
		if (VertexDegreeList[i] % 2 != 0) {
			check = false;
			break;
		}
	}
	return check;
}

bool UnorientedGraph::CheckHamiltonian()
{
	if (Nvertex == 2) {
		return false;
	}
	int mn = INT32_MAX;
	for (int i = 0; i < VertexDegreeList.size(); i++) {
		mn = min(mn, VertexDegreeList[i]);
	}
	return mn >= (Nvertex / 2);
}

void UnorientedGraph::MakeEulerian() {
	vector<int> Path;
	vector<int> odd;
	for (int i = 0; i < VertexDegreeList.size(); i++) {
		if (VertexDegreeList[i] % 2 != 0) {
			odd.push_back(i);
		}
	}
	cout << "--Deleted edges--" << endl;
	for (int i = 0; i < odd.size(); i++) {
		for (int j = i+1; j < odd.size(); j++) {
			if (AdjMatrix[i][j]) {
				AdjMatrix[i][j] = 0;
				AdjMatrix[j][i] = 0;
				VertexDegreeList[i]--;
				VertexDegreeList[j]--;

				cout << "(" << i << ", " << j << ")" << endl;
				break;
			}
		}
	}
	cout << endl;

	for (int i = 0; i < VertexDegreeList.size(); i++) {
		if (VertexDegreeList[i] % 2 != 0) {
			for (int j = i + 1; j < VertexDegreeList.size(); j++) {
				bool flag = false;
				if (VertexDegreeList[j] % 2 != 0) {
					for (int k = 0; k < VertexDegreeList.size(); k++) {
						if (AdjMatrix[i][k] and AdjMatrix[j][k] and VertexDegreeList[k] > 3) {
							AdjMatrix[i][k] = 0;
							AdjMatrix[k][i] = 0;

							AdjMatrix[j][k] = 0;
							AdjMatrix[k][j] = 0;

							VertexDegreeList[i]--;
							VertexDegreeList[j]--;
							VertexDegreeList[k] -= 2;
							flag = true;
							break;

						}
					}
					if (flag) {
						break;
					}
				}
			}
		}

	}
	
	AdjList = vector<vector<int>>(Nvertex, vector<int>());
	Nedge = 0;
	for (int i = 0; i < Nvertex; i++) {
		Nedge += VertexDegreeList[i];
		for (int j = 0; j < Nvertex; j++) {
			if (AdjMatrix[i][j]) {
				AdjList[i].push_back(j);
			}
			else {
				WeightMatrix[i][j] = INFINITY;
			}
		}
	}
	Nedge = Nedge / 2;
}

vector<int> UnorientedGraph::Fleury() {
	vector<int> res;
	stack<int> q;
	int v = 0;
	q.push(v);
	while (!q.empty())
	{
		v = q.top();
		int sumi = 0;
		for (int i = 0; i < Nvertex; i++)
		{
			sumi += AdjMatrix[v][i];
			if (sumi != 0)
			{
				break;
			}
		}
		if (sumi == 0)
		{
			res.push_back(v);
			q.pop();
		}
		else
		{
			for (int i = 0; i < Nvertex; i++)
			{
				if (AdjMatrix[v][i])
				{
					q.push(i);
					AdjMatrix[v][i] = 0;
					AdjMatrix[i][v] = 0;
					break;
				}
			}
		}
	}

	return res;
}

void UnorientedGraph::MakeHamiltonian() {
	// Check on your graph
	
	//AdjMatrix = InputMatrix<int>();
	//Nvertex = AdjMatrix.size();
	//AdjList = vector<vector<int>>(Nvertex, vector<int>());
	//VertexDegreeList = vector<int>(Nvertex, 0);
	//for (int i = 0; i < Nvertex; i++) {
	//	for (int j = 0; j < Nvertex; j++) {
	//		VertexDegreeList[i] += AdjMatrix[i][j];
	//		if (AdjMatrix[i][j]) {
	//			AdjList[i].push_back(j);
	//		}

	//	}
	//}

	int degree = Nvertex / 2;
	vector<int> vertices;
	for (int i = 0; i < VertexDegreeList.size(); i++) {
		if (VertexDegreeList[i] < degree) {
			vertices.push_back(i);
		}
	}
	for (const int v : vertices) {
		for (int i = 0; i < Nvertex; i++) {
			if (VertexDegreeList[v] >= degree) {
				break;
			}
			if (!AdjMatrix[v][i] and v != i) {
				AdjMatrix[v][i] = 1;
				AdjMatrix[i][v] = 1;
				AdjList[v].push_back(i);
				AdjList[i].push_back(v);
				WeightMatrix[v][i] = PascalsDistribution(7, 0.3, rand());
				WeightMatrix[i][v] = WeightMatrix[v][i];
				
				VertexDegreeList[v]++;
				VertexDegreeList[i]++;
			}
		}
	}
}

void UnorientedGraph::HamiltonianCycles(ofstream& out) {
	// Initially value of boolean
	// flag is false
	hasCycle = false;

	// Store the resultant path
	vector<int> path;
	path.push_back(0);

	for (int i = 0; i < Nvertex; i++)
		visited[i] = false;

	visited[0] = true;

	// Function call to find all
	// hamiltonian cycles
	double ans = INT32_MAX;
	FindHamCycle(AdjMatrix, 1, path, out,ans);

	if (!hasCycle) {

		// If no Hamiltonian Cycle
		// is possible for the
		// given graph
		cout << "No Hamiltonian Cycle" << "possible " << endl;
		return;
	}
	cout << "TSP answer: " << ans << endl << endl;
}
// Function to check if a vertex v
// can be added at index pos in
// the Hamiltonian Cycle
bool isSafe(int v, vector<vector<int>>& graph, vector<int>& path, int pos)
{

	// If the vertex is adjacent to
	// the vertex of the previously
	// added vertex
	if (graph[path[pos - 1]][v] == 0)
		return false;

	// If the vertex has already
	// been included in the path
	for (int i = 0; i < pos; i++)
		if (path[i] == v)
			return false;

	// Both the above conditions are
	// not true, return true
	return true;
}

// Recursive function to find all
// hamiltonian cycles
void UnorientedGraph::FindHamCycle(vector<vector<int>> graph, int pos, vector<int> path, ofstream& out, double& ans)
{
	// If all vertices are included
	// in Hamiltonian Cycle
	if (pos == Nvertex) {

		// If there is an edge
		// from the last vertex to
		// the source vertex
		if (graph[path[path.size() - 1]][path[0]] != 0) {

			// Include source vertex
			// into the path and
			// print the path
			path.push_back(0);

			 // TODO: заново пересоздаёт файл, надо в одном месте открыть его и всё
			double w = 0;
			for (int i = 0; i < path.size() - 1; i++) {
				w += WeightMatrix[path[i]][path[i+1]];
			}
			for (int i = 0; i < path.size(); i++) {
				cout << path[i] << " ";
				out << path[i] << " ";
			}
			if (w < ans) {
				ans = w;
			}
			cout << ", weight = " << w << endl;
			out << ", weight = " << w << endl;
			// Remove the source
			// vertex added
			path.pop_back();

			// Update the hasCycle
			// as true
			hasCycle = true;
		}
		return;
	}

	// Try different vertices
	// as the next vertex
	for (int v = 0; v < Nvertex; v++) {

		// Check if this vertex can
		// be added to Cycle
		if (isSafe(v, graph, path, pos) && !visited[v]) {

			path.push_back(v);
			visited[v] = true;

			// Recur to construct
			// rest of the path
			FindHamCycle(graph, pos + 1, path, out,ans);

			// Remove current vertex
			// from path and process
			// other vertices
			visited[v] = false;
			path.pop_back();
		}
	}
}