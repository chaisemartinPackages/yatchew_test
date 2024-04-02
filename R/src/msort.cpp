#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <Rcpp.h>

using namespace Rcpp;
using namespace std;

float norm(NumericVector V1, NumericVector V2) {
    float res = 0.0;
    for (int h = 0; h < V1.length(); h ++) {
        res += pow(V1[h] - V2[h], 2);
    }
    return(pow(res, 0.5));
}

bool notin(IntegerVector vec, int num) {
    for (int j = 0; j < vec.size(); j++) {
        if (vec[j] == num) {
            return(false);
        }
    }
    return(true);
}

int min_e_dist(NumericMatrix mat, int i, IntegerVector vec) {
    float dist = 0; int res = 0; float temp;
    for (int j = 0; j < mat.nrow(); j++) {
        temp = (i < j)? mat(j, i) : mat(i, j);
        if (((dist != 0) & (temp < dist)) | (dist == 0)) {
            if (notin(vec, j)) {
                dist = temp; res = j;
            }
        }
    }
    return(res);
}

// [[Rcpp::export]]
IntegerVector msort(NumericMatrix base) {
    const int nrow = base.nrow();

    // Preallocation of distances //
    NumericMatrix mat(nrow, nrow);
    for (int i = 1; i < nrow; i ++) {
        for (int j = 0; j < i; j++) {
            mat(i, j) = norm(base(i,_), base(j,_));
        }
    }

    IntegerVector id(nrow), checked;
    int start = 0;
    id[0] = start; 
    int S = start;

    for (int i = 1; i < nrow; i++) {
        checked.push_back(S);
        S = min_e_dist(mat, S, checked);
        id[i] = S;
    }
    
    for (int j = 0; j < nrow; j++) {
        id[j] += 1;
    }
    return(id);
}