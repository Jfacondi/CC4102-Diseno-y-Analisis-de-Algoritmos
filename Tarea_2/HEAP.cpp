#include <iostream>
#include <vector>
#include <queue>
#include <limits>
#include <utility>
#include <random>
#include <set>
#include <ctime>
#include <cstdlib>
#include <cmath> // Para la función pow

using namespace std;

// Estructura que representa un nodo en el grafo
struct Node {
    long long vertex;  // Vértice del nodo
    double distance;   // Distancia del nodo desde el origen
    long long previous; // Nodo previo en el camino mínimo

    // Constructor del nodo
    Node(long long v, double d, long long p) : vertex(v), distance(d), previous(p) {}
};

// Clase que implementa un MinHeap (montículo mínimo)
class MinHeap {
private:
    vector<Node> heap;      // Vector que almacena los nodos en el montículo
    vector<int> position;   // Vector que almacena las posiciones de los nodos en el montículo

    // Función para mantener la propiedad del montículo mínimo hacia abajo
    void heapifyDown(int idx) {
        int smallest = idx;
        int left = 2 * idx + 1;
        int right = 2 * idx + 2;

        if (left < heap.size() && heap[left].distance < heap[smallest].distance)
            smallest = left;

        if (right < heap.size() && heap[right].distance < heap[smallest].distance)
            smallest = right;

        if (smallest != idx) {
            swap(position[heap[idx].vertex], position[heap[smallest].vertex]);
            swap(heap[idx], heap[smallest]);
            heapifyDown(smallest);
        }
    }

    // Función para mantener la propiedad del montículo mínimo hacia arriba
    void heapifyUp(int idx) {
        if (idx && heap[(idx - 1) / 2].distance > heap[idx].distance) {
            swap(position[heap[idx].vertex], position[heap[(idx - 1) / 2].vertex]);
            swap(heap[idx], heap[(idx - 1) / 2]);
            heapifyUp((idx - 1) / 2);
        }
    }

public:
    // Constructor de MinHeap
    MinHeap(int size) {
        heap.reserve(size);    // Reserva espacio en el vector de heap
        position.resize(size, -1);  // Inicializa el vector de posiciones
    }

    // Verifica si el montículo está vacío
    bool isEmpty() const {
        return heap.empty();
    }

    // Inserta un nodo en el montículo
    void insert(Node node) {
        heap.push_back(node);
        int idx = heap.size() - 1;
        position[node.vertex] = idx;
        heapifyUp(idx);
    }

    // Extrae el nodo con la distancia mínima del montículo
    Node extractMin() {
        if (isEmpty()) throw runtime_error("Heap underflow");
        Node root = heap[0];    // Nodo raíz (mínimo)
        Node lastNode = heap.back();   // Último nodo en el montículo
        heap[0] = lastNode;
        position[lastNode.vertex] = 0;
        heap.pop_back();
        position[root.vertex] = -1;
        heapifyDown(0);
        return root;
    }

    // Disminuye la distancia de un nodo en el montículo
    void decreaseKey(int vertex, double newDist) {
        int idx = position[vertex];
        heap[idx].distance = newDist;
        heapifyUp(idx);
    }

    // Verifica si el montículo contiene un vértice específico
    bool contains(int vertex) const {
        return position[vertex] != -1;
    }
};
