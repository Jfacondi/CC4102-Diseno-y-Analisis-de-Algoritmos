#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <algorithm>

using namespace std;



//B es el tamaño de un bloque

#include <set>
#include <algorithm>

// pointer to 0xbaadf00dbaadf00d


struct Point{
    double x;
    double y;
    //constructor
    Point(double x, double y) : x(x), y(y) {}

    // Define operator== for Point class
    bool operator==(const Point& other) const {
        return x == other.x && y == other.y;
    }
};

// Forward declaration
class MTreeNode;

// Entry
class Entry {
public:
    Point point;
    double coveringRadius;
    MTreeNode* child; // a
    //setters
    void setchild(MTreeNode* a){
        child = a;
    }
    // Default constructor with default child
    Entry(Point p, double cr, MTreeNode* a) : point(p), coveringRadius(cr), child(a) {}
    // Default constructor with default child
    Entry(Point p) : point(p), coveringRadius(0), child(nullptr) {}

};

int B = 4096/sizeof(Entry); // Tamaño de un bloque de disco en bytes
int b = B/2; //tamaño minimo de un nodo
int max_node_size = B; //tamaño maximo de un nodo

// M-tree Node
class MTreeNode {
public:

    //bool isFather;
   
    // Add entry
    vector<Entry> entries;
    bool isLeaf = false;
    // Constructor
    MTreeNode() {};
    // Constructor who receives a node
    MTreeNode(MTreeNode* node) {
        entries = node->entries;
    }
    void addEntry(Point point, double cr, MTreeNode* childPage) {
        entries.push_back(Entry(point, cr, childPage));
    }
    // Add entry with default child
    void addEntry(Point point) {
        entries.push_back(Entry(point));
    }

};

// Distancia euclidiana entre dos puntos
double euclidean_dist(Point p1, Point p2) {
    return sqrt(pow(p1.x - p2.x, 2) + pow(p1.y - p2.y, 2));
}
//calcula la altura de un arbol
int height(MTreeNode* node, int h=0) {
    if (node == nullptr) {
        return 0;
    }
    if (node->isLeaf ) {
        return 1;
    }
    else {
         for(Entry &e : node->entries){
            int h1 = height(e.child, h);
            h = max(h, h1);
            return h + 1;
    }
    }
    return h;
}


void rellenarT_Sup(MTreeNode* T_Sup, vector<MTreeNode*> T_prime, vector<Point> samples){
    if(T_Sup == nullptr){
        return;
    }
    if(T_Sup->isLeaf){
        for(Entry &e : T_Sup->entries){
            for(int i = 0; i < samples.size(); i++){
                if(e.point == samples[i]){
                    e.child = T_prime[i];
                    break; 
                }
            }
        }
    }
    else{
        for(Entry &e : T_Sup->entries){
            rellenarT_Sup(e.child, T_prime, samples);
        }
    }
}


void radiosCobertores(MTreeNode* T){
    for(Entry &e : T->entries){
        if(e.child==nullptr || e.child->entries.empty()){
            continue;
        }
        double maxDist = 0;
        for(Entry &e1 : e.child->entries){
            maxDist = max(maxDist, euclidean_dist(e.point, e1.point));
        }
        e.coveringRadius = maxDist;
        radiosCobertores(e.child);
    }
}

void mtreeheightfinder(MTreeNode* T, int h, vector<MTreeNode*>& T_find, vector<Point>& samples){ 
    if(T == nullptr){
        return;
    }
    if(T->isLeaf){
        if(height(T) == h){
            T_find.push_back(T);
            for(Entry &e : T->entries){
                samples.push_back(e.point);
                break;
            }
        }
    }
    else{
        for(Entry &e : T->entries){
            if (height(e.child) == h){
                T_find.push_back(e.child);
                samples.push_back(e.point);
            }
            else if (height(e.child) < h){
                mtreeheightfinder(e.child, h, T_find, samples);
            }
        }
    }
}
void mtreebooleansetter(MTreeNode* T){
    if(T == nullptr){
        return;
    }
    else if (T->entries.empty() || T->entries[0].child == nullptr){
        T->isLeaf = true;
    }
    else{
        T->isLeaf = false;
        for(Entry &e : T->entries){
            e.child->isLeaf = false;
            mtreebooleansetter(e.child);
        }
    }
}
MTreeNode* CPMethod(vector<Point>& P) {
    // Step 1
    if (P.size() <= B) {
        // Se crea un árbol T, se insertan todos los puntos a T y se retorna T
        MTreeNode* T = new MTreeNode(); // Nodo hoja
        for (int i = 0; i < P.size(); i++) {
            T->addEntry(P[i]);
            T->isLeaf = true;
        }
        return T;
    }
    // Step 2
    //De manera aleatoria se eligen k = min(B, n/B) puntos de P, que los llamaremos samples pf1, . . . , pfk
    // takes the ceil of n/B with n being the size of P
    int b = ceil((double)P.size() / B);
    int k = min(B, b);
    vector<Point> samples;
    set<int> set_idxs;
    vector<vector<Point> > F(k);
        for (int i = 0; i < k; i++) {
        // Se eligen k = min(B, n/B) puntos de P, que los llamaremos samples pf1, . . . , pfk
        int idx = rand() % P.size();
        while (set_idxs.find(idx) != set_idxs.end()) {
            // Si ya se eligió ese punto, se elige otro
            idx = rand() % P.size();
        }
        samples.push_back(P[idx]); //F
        set_idxs.insert(idx);        
    }
    // Step 3
    //Se le asigna a cada punto en P un punto en sample, que es el más cercano. Con eso se puede construir k conjuntos F1, . . . , Fk. // F1, F2, ..., Fk
    for (int i = 0; i < P.size(); i++) {
        double minDist = 1e9;
        int minIdx = -1;
        for (int j = 0; j < k; j++) {
            double dist = euclidean_dist(P[i], samples[j]);
            if (dist < minDist) {
                minDist = dist;
                minIdx = j;
            }
        }
        F[minIdx].push_back(P[i]);
    }
    //step 4 Si algún Fj es tal que |Fj| < b, quitamos pfj de samples y trabajamos con los puntos de Fj para redistribuirlos en los conjuntos F1, . . . , Fk.
    for (int i = 0; i < F.size(); i++) {
        if (F[i].size() < b) {
            
            // Eliminar el punto pfj de samples
            // Trabajar con los puntos de Fj para redistribuirlos en los conjuntos F1, . . . , Fk
            for (int j = 0; j < F[i].size(); j++) {
                double minDist = 1e9;
                int minIdx = -1;
                for (int l = 0; l < samples.size(); l++) {
                    if (l == i) {
                        continue;
                    }
                    else{
                    double dist = euclidean_dist(F[i][j], samples[l]);
                    if (dist < minDist) {
                        minDist = dist;
                        minIdx = l;
                    }
                    }
                }
                F[minIdx].push_back(F[i][j]);
            }
            samples.erase(samples.begin() + i);
            F.erase(F.begin() + i);
        }
    }
    //Step 5
    if (samples.size() == 1) {
        return CPMethod(P);
    }
    // Step 6
    //Se realiza recursivamente el algoritmo CP en cada Fj, obteniendo el árbol Tj
    vector<MTreeNode*> T;
    for (int i = 0; i < F.size(); i++) {
        MTreeNode* node = new MTreeNode(CPMethod(F[i]));
        T.push_back(node);
    }
    // Step 7
    for (int i = 0; i < T.size(); i++) {
        if (T[i]->entries.size() < b) {
            //cout << "Tamaño menor a b es \n" << T[i]->entries.size() << endl;
            // Eliminar el punto pfj de F
            samples.erase(samples.begin() + i);
            // Agregar los puntos pertinentes de T[i] a F
            for (int j = 0; j < T[i]->entries.size(); j++) {
                if (T[i]->entries[j].child == nullptr) {
                    continue;
                }
                samples.push_back(T[i]->entries[j].point);
                T.push_back(T[i]->entries[j].child);
            }
            // Eliminar T[i] y agregar sus hijos a T
            T.erase(T.begin() + i);
        }
    }
    // Step 8
    //Etapa de balanceamiento: Se define h como la altura mínima de los árboles Tj
    //. Se define T' inicialmente como un conjunto vacío.
    
    int h = 1e9;
    for (int i = 0; i < T.size(); i++) {
        h = min(h, height(T[i]));
    }
    vector<MTreeNode*> T_prime;

    for (int i = 0; i < T.size(); i++) {
        if (height(T[i]) == h) {
            //put in T' in index i
            T_prime.insert(T_prime.begin() + i, T[i]);
        }
        else {
            samples.erase(samples.begin() + i);
            vector<MTreeNode*> T_find;
            vector<Point> samples2;
            mtreeheightfinder(T[i], h, T_find, samples2);
            for (int j = 0; j < T_find.size(); j++) {
                T_prime.push_back(T_find[j]);    
                samples.push_back(samples2[j]);
            }
        }
    }
    // Step 10
    MTreeNode* T_Sup = CPMethod(samples);
    // Step 11
    //Se añade T' a TSup
    rellenarT_Sup(T_Sup, T_prime, samples);
    mtreebooleansetter(T_Sup);
    radiosCobertores(T_Sup);
    // Step 13
    return T_Sup;

}


void mtreeprinter(MTreeNode* T){
    if(T == nullptr){
        return;
    }
    if(T->isLeaf){
        for(Entry &e : T->entries){
            cout << "end " << e.point.x << " " << e.point.y << endl;
        }
    }
    else{
        for(Entry &e : T->entries){
            cout << "entry size " << e.child->entries.size() << endl;
            mtreeprinter(e.child);
        }
    }
}

Point medoide(vector<Point>& cluster){
    double minDist = 1e9;
    Point medoide = Point(0, 0);
    for (int i = 0; i < cluster.size(); i++) {
        double sum = 0;
        for (int j = 0; j < cluster.size(); j++) {
            sum += euclidean_dist(cluster[i], cluster[j]);
        }
        if (sum < minDist) {
            minDist = sum;
            medoide = cluster[i];
        }
    }
    return medoide;
}


//funcion para sacar el par mas cercano de clusters

//el par mas cercano es un par de clusters c1, c2 tal que dist(c1, c2) ≤ dist(ci, cj ) para cualquier otro par ci, cj .
//input: un conjunto de clusters C
//output: un par de clusters c1, c2
pair<vector<Point>, vector<Point>> parMasCCercano(vector<vector<Point>>& C) {
    double minDist = 1e9;
    pair<vector<Point>, vector<Point>> par;
    //optimización: no es necesario comparar todos los pares de clusters
    vector<Point> medoides;
    for(vector<Point> optpunt : C){
        medoides.push_back(medoide(optpunt));
    }

    for (int i = 0; i < medoides.size(); i++) {
        for (int j = i + 1; j < medoides.size(); j++) {
            double dist = euclidean_dist(medoides[i], medoides[j]);
            if (dist < minDist) {
                minDist = dist;
                // Envolver los vectores en un par antes de asignarlos
                par = make_pair(C[i], C[j]);
            }
        }
    }
    return par;
}

//funcion para sacar el vecino mas cercano de un cluster:
//el vecino más cercano es otro cluster c′ tal que no hay otro cluster que su distancia a c sea menor a la distancia entre c y c′

//input: un conjunto de clusters C, un cluster c
//output: el vecino más cercano de c

vector<Point> vecinoMasCercano(vector<vector<Point>>& C, vector<Point>& c) {
    double minDist = 1e9;
    vector<Point> vecino;

    //optimizacion:
    for (int i = 0; i < C.size(); i++) {
        if (C[i] == c) {
            continue;
        }
        double dist = euclidean_dist(medoide(c), medoide(C[i]));
        if (dist < minDist) {
            minDist = dist;
            vecino = C[i];
        }
    }
    return vecino;
}

//radio cobertor maximo entre dos grupos resultantes
// Función para calcular el radio cobertor máximo entre dos grupos de puntos dados los medoides
double calculateMaxCoveringRadius(vector<Point>& group1, vector<Point>& group2, Point centroid1, Point centroid2) {
    double maxRadius = 0.0;
    // Itera sobre todos los puntos de group1 y group2 para calcular la distancia máxima a sus centroides
    for (const auto& p : group1) {
        maxRadius = max(maxRadius, euclidean_dist(centroid1, p));
    }
    for (const auto& p : group2) {
        maxRadius = max(maxRadius, euclidean_dist(centroid2, p));
    }
    return maxRadius;
}

// Función para dividir un conjunto en dos grupos usando la política MinMax split
pair<vector<Point>, vector<Point>> splitMinMax(vector<Point>& points) {
    double minMaxRadius = numeric_limits<double>::max();
    pair<vector<Point>, vector<Point>> optimalSplit;
    // Itera sobre todos los pares de puntos en el conjunto
    for (size_t i = 0; i < points.size(); ++i) {
        for (size_t j = i + 1; j < points.size(); ++j) {
            // Calcula el radio cobertor máximo entre los grupos resultantes si se divide el conjunto en dos partes usando los puntos i y j como centroides
            double maxRadius = calculateMaxCoveringRadius(points, points, points[i], points[j]);
            // Actualiza el par de puntos que minimiza el radio cobertor máximo
            if (maxRadius < minMaxRadius) {
                minMaxRadius = maxRadius;
                optimalSplit.first.clear();
                optimalSplit.second.clear();
                // Divide el conjunto en dos grupos usando los puntos i y j como centroides
                for (size_t k = 0; k < points.size(); ++k) {
                    double distToCentroid1 = euclidean_dist(points[i], points[k]);
                    double distToCentroid2 = euclidean_dist(points[j], points[k]);
                    if (distToCentroid1 < distToCentroid2) {
                        optimalSplit.first.push_back(points[k]);
                    } else {
                        optimalSplit.second.push_back(points[k]);
                    }
                }
            }
        }
    }
    return optimalSplit;
}

//funcion Cluster, que retorna un conjunto de clusters de tamaño entre b y B.


vector<vector<Point>> Cluster(vector<Point>& Cin) {
    vector<vector<Point>> Cout;
    vector<vector<Point>> C;
    for (int i = 0; i < Cin.size(); i++) {
        C.push_back({Cin[i]});
    }
    //3: Mientras |C| > 1:
    while (C.size() > 1) {
        printf("size of C on Cluster:%d\n", C.size());

        pair<vector<Point>, vector<Point>> par = parMasCCercano(C);
        vector<Point> c1;
        vector<Point> c2;

        if (par.first.size() > par.second.size() ){
            c1 = par.first;
           c2 = par.second;
        } else{
          c1 = par.second;
          c2 = par.first;
        }
       
            if (c1.size() + c2.size() <= B) {
                C.erase(remove(C.begin(), C.end(), c1), C.end());
                C.erase(remove(C.begin(), C.end(), c2), C.end());
                //se añade c1 ∪ c2 a C
                vector<Point> c1uc2;
                c1uc2.clear();
                c1uc2.insert(c1uc2.end(), par.first.begin(), par.first.end());
                c1uc2.insert(c1uc2.end(), par.second.begin(), par.second.end());
                C.push_back(c1uc2); 
            }
        else{
            //se remueve c1 de C y se añade c1 a Cout
            C.erase(remove(C.begin(), C.end(), par.first), C.end());
            Cout.push_back(par.first);
        }
    
}
    //4: Sea c el último elemento de C
    vector<Point> c = C.back();
    //5: Si |Cout| > 0:
    vector<Point> c_prime;
    if (Cout.size() > 0) {
        //5.1 definimos c′ como el vecino más cercano a c en Cout. Removemos c′ de Cout
        vector<Point> c_prime = vecinoMasCercano(Cout, c);
        Cout.erase(remove(Cout.begin(), Cout.end(), c_prime), Cout.end());
        //5.2 Si no, se define c′ = {}.
    }
   
    //6: Si |c ∪ c′ | ≤ B:
    //armamos c_uc_prime esta en su prime
    vector<Point> c_uc_prime;
    c_uc_prime.clear();
    c_uc_prime.insert(c_uc_prime.end(), c.begin(), c.end());
    c_uc_prime.insert(c_uc_prime.end(), c_prime.begin(), c_prime.end());

    if(c.size() + c_prime.size() <= B){
        //6.1 Añadimos c ∪ c′ a Cout.
       
        Cout.push_back(c_uc_prime);

    }
    else{
        //6.2 Si no, dividimos c ∪ c′ en c1 y c2 usando MinMax split policy. Se añaden c1 y c2 a Cout.
        // MinMax split policy: Se considera todos los posibles pares de puntos, y alternadamente se van agregando el punto más cercano a alguno de estos centros
        // (esto garantiza que la división sea balanceada) y se calcula el radio cobertor máximo entre estos dos grupos resultantes.
        // Esto se prueba para todo par de puntos y se elige el par que tenga el mínimo radio cobertor máximo.
        pair<vector<Point>, vector<Point>> split = splitMinMax(c_uc_prime); 
        Cout.push_back(split.first);
        Cout.push_back(split.second);
    }
    return Cout;
}

//funcion OutputHoja
//input: un conjunto de puntos Cin. 
//output: una tupla (g, r, a) donde g es el medoide primario de Cin, r es llamado el radio cobertor y a la dirección del hijo respectivo.
/*
OutputHoja: Retorna una tupla (g, r, a) donde g es el medoide primario de Cin, r es llamado el
radio cobertor y a la dirección del hijo respectivo.
Input: Cin
1. Sea g el medoide primario de Cin. Sea r = 0. Sea C = {} (el que corresponderá al nodo
hoja).
2. Por cada p ∈ Cin: Añadimos (p, null, null) a C. Seteamos r = max(r, dist(g, p))
3. Guardamos el puntero a C como a
4. Retornamos (g, r, a)

*/



//2da version de OutputHoja

   Entry OutputHoja(vector<Point>& Cin){
    Point g = medoide(Cin);
    double r = 0;
    //C es un nodo hoja
    MTreeNode* C = new MTreeNode(); 
    C->isLeaf = true;
    for (int i = 0; i < Cin.size(); i++) {
        C->addEntry(Cin[i]);
        r = max(r, euclidean_dist(g, Cin[i]));
    }

    //3: Guardamos el puntero a C como a

    Entry e = Entry(g, r, C);
   
    //Retornamos (g, r, a)
    return e;
   }

/*
OutputInterno: Retorna (G, R, A) donde G es el medoide primario del conjunto de puntos
Cin = {g|∃(g, r, a) ∈ Cmra}, R el radio cobertor, y A la dirección del hijo respectivo.
input: Cmra, un conjunto de tuplas (g, r, a) retornadas por OutputHoja
*/

//input: un conjunto de tuplas Cmra (en verdad es un conjunto de entradas)
//output: una tupla (G, R, A) donde G es el medoide primario del conjunto de puntos Cin = {g|∃(g, r, a) ∈ Cmra}, R el radio cobertor, y A la dirección del hijo respectivo.

 Entry OutputInterno(vector<Entry>& Cmra){
    vector<Point> Cin;
    for (int i = 0; i < Cmra.size(); i++) {
        Cin.push_back(Cmra[i].point);
    }
    Point G = medoide(Cin);
    double R = 0;
    MTreeNode* C = new MTreeNode();
    C->isLeaf = false;
    for (int i = 0; i < Cmra.size(); i++) {
        C->addEntry(Cmra[i].point, Cmra[i].coveringRadius, Cmra[i].child);
        R = max(R, euclidean_dist(G, Cmra[i].point) + Cmra[i].coveringRadius);
    }
    //3: Guardamos el puntero a C como A.

    Entry e = Entry(G, R, C);
    //Retornamos (G, R, A)
    return e;
}


//funcion SS
//input: un conjunto de puntos Cin
//output: la raíz del M-tree construído.

MTreeNode SS(vector<Point>& Cin){
    //1: Si |Cin| ≤ B: Se define (g, r, a) = OutputHoja(Cin) y se retorna a.
    if (Cin.size() <= B) {
        Entry e = OutputHoja(Cin);
        //queremos retonar solo a.
        return e.child;
    }
    //2: Sea Cout = Cluster(Cin). Sea C = {}.
    vector<vector<Point>> Cout = Cluster(Cin);
    vector<Entry> C;
   
    //3: Por cada c ∈ Cout: Se añade OutputHoja(c) a C
    for (int i = 0; i < Cout.size(); i++) {
        C.push_back(OutputHoja(Cout[i]));
    }
    //4: Mientras |C| > B:
    while (C.size() > B) {
        printf("size of C %d\n", C.size());
        //4.1 Sea Cin = {g|(g, r, a) ∈ C}. Sea Cout = Cluster(Cin). Sea Cmra = {}.
        vector<Point> Cin;
        for (int i = 0; i < C.size(); i++) {
            Cin.push_back(C[i].point);
        }
        vector<vector<Point>> Cout = Cluster(Cin);
        vector<Entry> Cmra;
        //4.2 Por cada c ∈ Cout: Sea s = {(g, r, a)|(g, r, a) ∈ C ∧ g ∈ c}, se añade s a Cmra
        for (int i = 0; i < Cout.size(); i++) {
            vector<Entry> s;
            //Cout es un conjunto de clusters. c es un cluster
            //C es un conjunto de entradas.
            //Para cada cluster c en Cout, se añade a s las entradas de C que tienen como punto a un punto en c
            for (int j = 0; j < C.size(); j++) {
                if (find(Cout[i].begin(), Cout[i].end(), C[j].point) != Cout[i].end()) {
                    s.push_back(C[j]);
                }
            }
            Cmra.insert(Cmra.end(), s.begin(), s.end());
        }

        
        //4.3 Sea C = {}.
        C.clear();
        //4.4 Por cada s ∈ Cmra: Añadir OutputInterno(s) a C
        for (int i = 0; i < Cmra.size(); i++) {
            C.push_back(OutputInterno(Cmra));
        }
    }
    
    //5: Sea (g, r, a) = OutputInterno(C)
    auto x = OutputInterno(C);
    //6: Se retorna a
    return x.child;

}

