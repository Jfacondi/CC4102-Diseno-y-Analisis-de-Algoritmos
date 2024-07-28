#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <algorithm>
#include <fstream>
#include <string>
#include "METHODS.cpp"
using namespace std;
//import the CPP_METHOD.cpp file
void search(MTreeNode* T, Point q, double r, vector<Point>& result, int& memory_access) {
    memory_access++;
    if (T->isLeaf) {
        for (Entry &e : T->entries) {
            if (euclidean_dist(e.point, q) <= r) {
                result.push_back(e.point);
            }
        }
    }
    else {
        for (Entry &e : T->entries) {
            if (euclidean_dist(e.point, q) <= r + e.coveringRadius) {
                if (euclidean_dist(e.point, q) <= r) {
                    result.push_back(e.point);
                }
                search(e.child, q, r, result, memory_access);
            }
        }
    }
}

int test(int i) {
    printf("Test %d\n", i);
    srand(time(0));
    vector<Point> P;
    int n = pow(2,i);
    for (int i = 0; i < n; i++) {
        double x = ((double) rand() / (RAND_MAX)) ;
        double y = ((double) rand() / (RAND_MAX)) ;
        P.push_back(Point(x, y));
    }
    clock_t start2, end2;
    start2 = clock();
    MTreeNode* T = CPMethod(P);
    end2 = clock();
    double time2 = (double)(end2 - start2) / CLOCKS_PER_SEC;
    vector<Point> X;
    for (int i = 0; i < 100; i++) {
        double x = ((double) rand() / (RAND_MAX)) ;
        double y = ((double) rand() / (RAND_MAX)) ;
        X.push_back(Point(x, y));
    }
    vector<int> vect_mem;
    double count = 0;
    int memory_access = 0;
    int memory_access_min = 1000000000;
    int memory_access_max = 0;
    int memoory = 0;
    double time = 0;
    clock_t start, end;

    for (Point &p : X) {
        vector<Point> Y;
        start = clock();
        search(T, p, 0.02, Y, memoory);
        end = clock();
        time += (double)(end - start) / CLOCKS_PER_SEC;
        count += (double)Y.size()*100/n;
        memory_access += memoory;
        vect_mem.push_back(memoory);
        if (memoory < memory_access_min) {
            memory_access_min = memoory;
        }
        if (memoory > memory_access_max) {
            memory_access_max = memoory;
        }
        memoory = 0;
    }
    count = count / X.size();
    time = time / X.size();
    memory_access = memory_access / X.size();
    int sum = 0;
    for(int i = 0; i < vect_mem.size(); i++) {
        sum = sum + (vect_mem[i] - memory_access) * (vect_mem[i] - memory_access);
    }
    double standard_deviation = sqrt(sum / vect_mem.size());
    double cot_inf = memory_access - 1.96 * standard_deviation / sqrt(vect_mem.size());
    double cot_sup = memory_access + 1.96 * standard_deviation / sqrt(vect_mem.size());
    //SAVE IN Resultados.txt the count, i, memory access and the time of search    
    FILE *f = fopen("Resultados.txt", "a");
    fprintf(f, "CPmethod : points percent %f, n pow %d, memory access %d , time of search %f, min memory access: %d, max memory access %d, cot inf %f, cot sup %f, execution time %f\n", count, i, memory_access, time, memory_access_min, memory_access_max, cot_inf, cot_sup, time2);
    if (i < 10) {
        start2 = clock();
        MTreeNode T2 = SS(P);
        end2 = clock();
        double time3 = (double)(end2 - start2) / CLOCKS_PER_SEC;
        vector<int> vect_mem2;
        double count2 = 0;
        int memory_access2 = 0;
        int memory_access_min2 = 1000000000;
        int memory_access_max2 = 0;
        int memoory2 = 0;
        double time2 = 0;
        for (Point &p : X) {
            vector<Point> Y;
            start = clock();
            search(&T2, p, 0.02, Y, memoory2);
            end = clock();
            time2 += (double)(end - start) / CLOCKS_PER_SEC;
            count2 += (double)Y.size()*100/n;
            memory_access2 += memoory2;
            vect_mem2.push_back(memoory2);
            if (memoory2 < memory_access_min2) {
                memory_access_min2 = memoory2;
            }
            if (memoory2 > memory_access_max2) {
                memory_access_max2 = memoory2;
            }
            memoory2 = 0;
        }
        count2 = count2 / X.size();
        time2 = time2 / X.size();
        memory_access2 = memory_access2 / X.size();
        int sum2 = 0;
        for(int i = 0; i < vect_mem2.size(); i++) {
            sum2 = sum2 + (vect_mem2[i] - memory_access2) * (vect_mem2[i] - memory_access2);
        }
        double standard_deviation2 = sqrt(sum2 / vect_mem2.size());
        double cot_inf2 = memory_access2 - 1.96 * standard_deviation2 / sqrt(vect_mem2.size());
        double cot_sup2 = memory_access2 + 1.96 * standard_deviation2 / sqrt(vect_mem2.size());
        fprintf(f, "SSmethod : points percent %f, n pow %d, memory access %d , time of search %f, min memory access: %d, max memory access %d, cot inf %f, cot sup %f, exec time %f\n", count2, i, memory_access2, time2, memory_access_min2, memory_access_max2, cot_inf2, cot_sup2, time3);
    }
    fclose(f);
    return 0;
}
int main() {
    for (int i = 10; i < 26; i++) {
        test(i);
    }
    return 0;
}
