#include <bits/stdc++.h>
#include "fibo.cpp"
#include "HEAP.cpp"

tuple<vector<double>, vector<long long>> DijkstraHeap(long long V, vector<vector<pair<long long, double>>> &adj, long long source) {
    vector<double> dist(V, numeric_limits<double>::infinity());
    vector<long long> prev(V, -1);
    dist[source] = 0;
    Node node = {source, 0, -1};
    MinHeap minHeap(V);
    minHeap.insert(node);

    while (!minHeap.isEmpty()) {
        Node u = minHeap.extractMin();
        for (auto &v : adj[u.vertex]) {
            if (dist[u.vertex] + v.second < dist[v.first]) {
                dist[v.first] = dist[u.vertex] + v.second;
                prev[v.first] = u.vertex;
                Node node = {v.first, dist[v.first], u.vertex};
                minHeap.insert(node);
            }
        }
    }
    return {dist, prev};
}

tuple <vector<double>,vector<long long>> DijkstraFibo(long long V, vector<vector<pair<long long, double>>> &adj, long long source) {
    vector<double> dist(V);
    vector<long long> prev(V);
    long long prevVertex = source;
    prev[source] = prevVertex;
    dist[source] = 0;
    FibonacciHeap fiboHeap(V);
    fiboHeap.push(source, 0);
    for (long long i = 0; i < V; i++) {
        if (i == source) continue;
        fiboHeap.push(i, numeric_limits<double>::infinity());
    }
    while (fiboHeap.totalNodes > 0) {
        NodeInfo u = fiboHeap.pop();
        long long uVertex = u.neighbor;
        for (auto &v : adj[uVertex]) {
            long long vVertex = v.first;
            double weight = v.second;
            if (dist[uVertex] + weight < dist[vVertex]) {
                dist[vVertex] = dist[uVertex] + weight;
                prev[vVertex] = uVertex;
                fiboHeap.push(vVertex, dist[vVertex]);
            }
        }
        


    }
    return {dist, prev};
}

void build_graph(vector<vector<pair<long long, double>>> &G, int V, int E){
    for (long long k = 0; k < V; k++) {
            long long rand_neigh = rand() % (k + 1); // ensure k is not 0
            double rand_weight = ((double) rand() / (RAND_MAX));
            G[k].push_back({rand_neigh, rand_weight});
            G[rand_neigh].push_back({k, rand_weight});
        }
        for (long long k = V; k < E; k++) {
            long long rand_neigh = rand() % V;
            long long rand_neigh2 = rand() % V;
            double rand_weight = ((double) rand() / (RAND_MAX));
            G[rand_neigh].push_back({rand_neigh2, rand_weight});
            G[rand_neigh2].push_back({rand_neigh, rand_weight});
        }
}

void test(int i, int j, int n_rep = 50) {
    long long v = pow(2, i);
    long long e = pow(2, j);
    clock_t t1, t2;
    double t_heap = 0;
    double t_fibo = 0;
    srand(time(0));
    for (int test_num = 1; test_num <= n_rep; test_num++) {
        printf("\rTest (i=%d, j=%d) %02d / %d", i, j, test_num, n_rep);
        cout.flush();
        vector<vector<pair<long long, double>>> adj(v);
        build_graph(adj, v, e);
        long long source = rand() % v;

        t1 = clock();
        auto [dist1, prev1] = DijkstraHeap(v, adj, source);
        t2 = clock();
        t_heap += (double)(t2 - t1) / CLOCKS_PER_SEC;

        t1 = clock();
        auto [dist2, prev2] = DijkstraFibo(v, adj, source);
        t2 = clock();
        t_fibo += (double) (t2 - t1) / CLOCKS_PER_SEC;
    }
    printf("\n");
    double mean_heap = t_heap / n_rep;
    double mean_fibo = t_fibo / n_rep;
    FILE *F = fopen("results.csv", "a");
    fprintf(F, "Heap,%d,%d,%f\n", i, j, mean_heap);
    fprintf(F, "Fibo,%d,%d,%f\n", i, j, mean_fibo);
    fclose(F);
}


int main() {
    FILE *F = fopen("results.csv", "w");
    fprintf(F, "queue,i,j,mean_t\n");
    fclose(F);
    for (int j = 16; j <= 22; j++) {
        test(10, j);
    }
    for (int j= 16; j <= 22; j++) {
        test(12, j);
    }
    for (int j = 16; j <= 22; j++) {
        test(14, j);
    }
    return 0;
}
