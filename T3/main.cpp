#include <iostream>
#include <string>
#include <functional>
#include <fstream>
#include <vector>
#include <unordered_set>
#include <random>
#include <chrono>
#include <numeric>
#include <algorithm>
#include <cmath>

using namespace std;

// Vectores para almacenar nombres de bebés y nombres de películas
vector<string> baby_names;
vector<string> film_names;

// Clase para una función de hash personalizada
class HashFunction {
public:
    // Constructor que inicializa el valor m
    HashFunction(size_t m) : m(m) {}
    
    // Sobrecarga del operador de función para calcular el valor hash
    size_t operator()(const string& key) const {
        hash<string> hash_fn;
        return hash_fn(key) % m;
    }
private:
    size_t m;  // Número de bits del Bloom Filter
};

// Clase para implementar un Bloom Filter
class BloomFilter {
public:
    size_t m, k;  // Número de bits del filtro y número de funciones hash
    vector<int> M;  // Vector de bits del filtro
    vector<HashFunction> hash_functions;  // Vector de funciones hash

    // Constructor por defecto
    BloomFilter() : m(0), k(0) {}

    // Constructor que inicializa el Bloom Filter
    BloomFilter(size_t m, size_t k, const vector<string>& names) : m(m), k(k) {
        M.resize(m, 0);  // Inicializa el vector de bits a 0
        for (size_t i = 0; i < k; i++) {
            hash_functions.push_back(HashFunction(m));  // Añade funciones hash
        }
        for (const string& name : names) {
            for (const HashFunction& h : hash_functions) {
                size_t index = h(name);  // Calcula el índice hash
                M[index] = 1;  // Marca el bit correspondiente en el filtro
            }
        }
    }

    // Método para buscar un nombre en el Bloom Filter
    bool search(const string& y) const {
        for (const HashFunction& h : hash_functions) {
            size_t index = h(y);  // Calcula el índice hash
            if (M[index] == 0) {
                return false;  // Si algún bit es 0, el nombre no está en el filtro
            }
        }
        return true;  // Si todos los bits son 1, el nombre puede estar en el filtro
    }
};

// Función para buscar un nombre en el vector de nombres de bebés
bool search(const string& y) {
    for (const string& name : baby_names) {
        if (name == y) return true;
    }
    return false;
}

// Función para buscar nombres sin usar Bloom Filter
vector<bool> busquedaSinFiltro(const vector<string>& names) {
    vector<bool> result(names.size(), false);
    for (size_t i = 0; i < names.size(); i++) {
        result[i] = search(names[i]);
    }
    return result;
}

// Función para buscar nombres usando Bloom Filter
vector<bool> busquedaConFiltro(const vector<string>& names, const BloomFilter& bf) {
    vector<string> names_to_search;
    for (const string& name : names) {
        if (bf.search(name)) {
            names_to_search.push_back(name);
        }
    }
    return busquedaSinFiltro(names_to_search);
}

// Función para generar un vector de nombres mezclando nombres de bebés y nombres de películas
vector<string> gen_vector_N(size_t N, double p) {
    vector<string> result;
    unordered_set<string> used_names;
    int num_baby_names = static_cast<int>(N * p);
    int num_film_names = N - num_baby_names;

    random_device rd;
    mt19937 g(rd());
    shuffle(baby_names.begin(), baby_names.end(), g);
    shuffle(film_names.begin(), film_names.end(), g);

    for (int i = 0; i < num_baby_names && i < baby_names.size(); ++i) {
        result.push_back(baby_names[i]);
        used_names.insert(baby_names[i]);
    }

    for (int i = 0; i < film_names.size() && result.size() < N; ++i) {
        if (used_names.find(film_names[i]) == used_names.end()) {
            result.push_back(film_names[i]);
            used_names.insert(film_names[i]);
        }
    }

    shuffle(result.begin(), result.end(), g);

    return result;
}

// Función para cargar nombres desde un archivo CSV
void load_csv(const string& filename, vector<string>& names) {
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Error al abrir el archivo\n";
        exit(1);
    }
    string line;
    getline(file, line); // Leer la cabecera si la tiene
    while (getline(file, line)) {
        if (!line.empty()) {
            names.push_back(line);
        }
    }
    file.close();
}

// Función de prueba para evaluar el rendimiento del Bloom Filter
void test(int i, double p) {
    ofstream f("results.txt", ios::app);
    f << "Test con N = 2^" << i << ", p = " << p << '\n';
    f.close();

    size_t N = pow(2, i);

    vector<double> tiempoSinFiltro, tiempoConFiltro, errores;

    for (int x = 1; x <= 50; x++) { // Cantidad de pruebas a realizar
        cout << "Test " << x << '\n';
        vector<string> names = gen_vector_N(N, p);

        auto start = chrono::high_resolution_clock::now();
        vector<bool> resultSinFiltro = busquedaSinFiltro(names);
        auto end = chrono::high_resolution_clock::now();
        chrono::duration<double> durationSinFiltro = end - start;
        tiempoSinFiltro.push_back(durationSinFiltro.count());

        size_t m = pow(2, 16); // Número de bits del Bloom Filter
        size_t k = 1000; // Número de funciones hash
        BloomFilter bf(m, k, baby_names);

        start = chrono::high_resolution_clock::now();
        vector<bool> resultConFiltro = busquedaConFiltro(names, bf);
        end = chrono::high_resolution_clock::now();
        chrono::duration<double> durationConFiltro = end - start;
        tiempoConFiltro.push_back(durationConFiltro.count());

        double error = 0;
        for (size_t j = 0; j < resultConFiltro.size(); j++) {
            if (j < resultSinFiltro.size() && resultConFiltro[j] != resultSinFiltro[j]) {
                error++;
            }
        }
        double pct_error = (error / resultConFiltro.size()) * 100;
        errores.push_back(pct_error);

       
    }

    double promedioTiempoSinFiltro = accumulate(tiempoSinFiltro.begin(), tiempoSinFiltro.end(), 0.0) / tiempoSinFiltro.size();
    double promedioTiempoConFiltro = accumulate(tiempoConFiltro.begin(), tiempoConFiltro.end(), 0.0) / tiempoConFiltro.size();
    double promedioError = accumulate(errores.begin(), errores.end(), 0.0) / errores.size();

    f.open("results.txt", ios::app);
    f << "Resultados promedio\n";
    f << "Caso numero:" << i << ' ' << p << '\n';
    f << "Tiempo sin filtro: " << promedioTiempoSinFiltro << '\n';
    f << "Tiempo con filtro: " << promedioTiempoConFiltro << '\n';
    f << "Porcentaje de error promedio: " << promedioError << "%\n";
    f.close();
}

// Función principal del programa
int main() {
    // Cargar los archivos CSV
    load_csv("Popular-Baby-Names-Final.csv", baby_names);
    load_csv("Film-Names.csv", film_names);

    // Ejecutar pruebas para diferentes valores de N y p
    for (int i = 10; i <= 16; i += 2) {
        for (double p = 0; p <= 1; p += 0.25) {
            test(i, p);
        }
    }
    return 0;
}
