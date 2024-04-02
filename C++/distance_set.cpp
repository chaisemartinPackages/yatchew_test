#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include "msort.h"

using namespace std;

float randn() {
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> dis(0.0, 1.0);
    return(dis(gen));
}


int main() {
    int nrow = 5;
    int ncol = 2;
    vector<vector<float>> vec(nrow, vector<float>(ncol));

    for (int i = 0; i < nrow; i++){
        for (int j = 0; j < ncol; j++) {
            vec[i][j] = randn();
        }
    }
    vector<int> ids;
    ids = msort(vec);
    print_vec(ids);
}