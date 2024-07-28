#include <bits/stdc++.h>

using namespace std;

// Clase que representa la información de un nodo vecino en el grafo
class NodeInfo {
public:
    long long neighbor; // Identificador del vecino
    double weight; // Peso del borde hacia el vecino
};

// Clase que representa un nodo en el Fibonacci Heap
class Node2 {
public:
    double key; // Clave del nodo (utilizada para comparar nodos)
    NodeInfo info; // Información del nodo
    int degree = 0; // Grado del nodo (número de hijos)
    bool mark = false; // Marca para operaciones de corte en cascada
    Node2 *child = nullptr; // Puntero al primer hijo
    Node2 *parent = nullptr; // Puntero al nodo padre
    Node2 *prev = nullptr; // Puntero al nodo anterior en la lista de hermanos
    Node2 *next = nullptr; // Puntero al nodo siguiente en la lista de hermanos
};

// Clase que representa un Fibonacci Heap
class FibonacciHeap {
public:
    vector<Node2*> nodes; // Vector que almacena los nodos
    Node2 *minNode = nullptr; // Puntero al nodo con la clave mínima
    int totalNodes = 0; // Total de nodos en el heap

    // Constructor que inicializa el heap con un tamaño dado
    FibonacciHeap(int size) {
        nodes = vector<Node2*>(size);
    }

    // Remueve un nodo de la lista de hermanos
    void removeNodeFromList(Node2 *node) {
        node->prev->next = node->next;
        node->next->prev = node->prev;
    }

    // Añade un nodo a la lista de hermanos
    void addNodeToList(Node2 *list, Node2 *nodeToAdd) {
        Node2 *temp = list->prev;
        nodeToAdd->next = list;
        nodeToAdd->prev = temp;
        list->prev = nodeToAdd;
        temp->next = nodeToAdd;
    }

    // Inserta un nuevo nodo en el heap
    void push(long long neighbor, double weight) {
        Node2 *newNode = new Node2();
        newNode->key = weight;
        newNode->info.neighbor = neighbor;
        newNode->info.weight = weight;
        nodes[neighbor] = newNode;
        totalNodes++;
        insertIntoHeap(newNode);
    }

    // Inserta un nodo en el heap y actualiza el minNode si es necesario
    void insertIntoHeap(Node2 *nodeToInsert) {
        if (minNode == nullptr) {
            nodeToInsert->prev = nodeToInsert;
            nodeToInsert->next = nodeToInsert;
            minNode = nodeToInsert;
        } else {
            addNodeToList(minNode, nodeToInsert);
            if (nodeToInsert->key < minNode->key) {
                minNode = nodeToInsert;
            }
        }
    }

    // Une dos Fibonacci Heaps
    void unionHeaps(FibonacciHeap *otherHeap) {
        if (otherHeap->minNode == nullptr) return;
        if (minNode == nullptr) {
            minNode = otherHeap->minNode;
            totalNodes = otherHeap->totalNodes;
            return;
        }

        // Combina las listas de raíces de ambos heaps
        Node2 *temp1 = minNode->next;
        Node2 *temp2 = otherHeap->minNode->prev;
        minNode->next = otherHeap->minNode;
        otherHeap->minNode->prev = minNode;
        temp2->next = temp1;
        temp1->prev = temp2;

        // Actualiza el nodo mínimo si es necesario
        if (otherHeap->minNode->key < minNode->key) {
            minNode = otherHeap->minNode;
        }
        totalNodes += otherHeap->totalNodes;
    }

    // Verifica si el heap está vacío
    bool isEmpty() {
        return minNode == nullptr;
    }

    // Extrae el nodo con la clave mínima del heap
    NodeInfo pop() {
        Node2 *min = minNode;
        if (min != nullptr) {
            // Añade los hijos del nodo mínimo a la lista de raíces
            if (min->child != nullptr) {
                Node2 *child = min->child;
                do {
                    Node2 *nextChild = child->next;
                    addNodeToList(minNode, child);
                    child->parent = nullptr;
                    child = nextChild;
                } while (child != min->child);
            }
            removeNodeFromList(min);
            if (min == min->next) {
                minNode = nullptr;
            } else {
                minNode = min->next;
                consolidate();
            }
            totalNodes--;
        }
        NodeInfo result = min->info;
        delete min;
        nodes[result.neighbor] = nullptr;
        return result;
    }

    // Consolidación del heap para mantener la estructura del Fibonacci Heap
    void consolidate() {
        double phi = (1 + sqrt(5)) / 2;
        int maxDegree = log(totalNodes) / log(phi);
        vector<Node2*> degreeTable(maxDegree + 1, nullptr);

        // Lista temporal de nodos raíz
        vector<Node2*> rootNodes;
        Node2 *current = minNode;
        do {
            rootNodes.push_back(current);
            current = current->next;
        } while (current != minNode);

        // Consolida nodos de la misma orden
        for (Node2 *node : rootNodes) {
            Node2 *x = node;
            int degree = x->degree;
            while (degreeTable[degree] != nullptr) {
                Node2 *y = degreeTable[degree];
                if (x->key > y->key) swap(x, y);
                linkNodes(y, x);
                degreeTable[degree] = nullptr;
                degree++;
            }
            degreeTable[degree] = x;
        }

        // Reconstruye la lista de raíces y actualiza el nodo mínimo
        minNode = nullptr;
        for (Node2 *node : degreeTable) {
            if (node != nullptr) {
                if (minNode == nullptr) {
                    minNode = node;
                    minNode->prev = minNode;
                    minNode->next = minNode;
                } else {
                    addNodeToList(minNode, node);
                    if (node->key < minNode->key) {
                        minNode = node;
                    }
                }
            }
        }
    }

    // Enlaza dos nodos en el heap
    void linkNodes(Node2 *y, Node2 *x) {
        removeNodeFromList(y);
        if (x->child == nullptr) {
            x->child = y;
            y->prev = y;
            y->next = y;
        } else {
            addNodeToList(x->child, y);
        }
        y->parent = x;
        y->mark = false;
        x->degree++;
    }

    // Disminuye la clave de un nodo específico
    void decreaseKey(long long neighbor, double newDistance) {
        Node2 *node = nodes[neighbor];
        if (node->key < newDistance) return;
        node->key = newDistance;
        node->info.weight = newDistance;
        Node2 *parent = node->parent;
        if (parent != nullptr && node->key < parent->key) {
            cutNode(node, parent);
            cascadingCut(parent);
        }
        if (node->key < minNode->key) {
            minNode = node;
        }
    }

    // Corta un nodo de su padre y lo añade a la lista de raíces
    void cutNode(Node2 *node, Node2 *parent) {
        parent->degree--;
        if (parent->degree == 0) {
            parent->child = nullptr;
        } else {
            if (parent->child == node) {
                parent->child = node->next;
            }
            removeNodeFromList(node);
        }
        addNodeToList(minNode, node);
        node->parent = nullptr;
        node->mark = false;
    }

    // Realiza el corte en cascada de un nodo
    void cascadingCut(Node2 *node) {
        Node2 *parent = node->parent;
        if (parent != nullptr) {
            if (!node->mark) {
                node->mark = true;
            } else {
                cutNode(node, parent);
                cascadingCut(parent);
            }
        }
    }
};