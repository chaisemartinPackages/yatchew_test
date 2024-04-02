#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <Rcpp.h>
using namespace Rcpp;

using namespace std;

// [[Rcpp::export]]
float norm(vector<float> V1, vector<float> V2) {
    float res = 0.0;
    for (int h = 0; h < V1.size(); h ++) {
        res += pow(V1[h] - V2[h], 2);
    }
    return(pow(res, 0.5));
}

// [[Rcpp::export]]
bool notin(vector<int> vec, int num) {
    for (int j = 0; j < vec.size(); j++) {
        if (vec[j] == num) {
            return(false);
        }
    }
    return(true);
}

// [[Rcpp::export]]
void print_mat(vector<vector<float>> mat, int xlim, int ylim) {
    for (int i = 0; i < xlim; i++) {
        for (int j = 0; j < ylim; j++) {
            cout << mat[i][j]  << " ";
        }
        cout << endl;
    }
}

// [[Rcpp::export]]
vector<int> msort(vector<vector<float>> vec) {
    const int nrow = vec.size();
    const int ncol = vec[0].size();

    // Preallocation of distances //
    vector<vector<float>> mat(nrow, vector<float>(nrow));
    for (int i = 1; i < nrow; i ++) {
        for (int j = 0; j < i; j++) {
            mat[i][j] = norm(vec[i], vec[j]);
        }
    }

    // Origin vector
    vector<float> origin(ncol);
    for (int j = 0; j < ncol; j ++) {
        origin[j] = 0.0;
    }

    // Preallocate the id set
    vector<int> id(nrow);
    int start = 0; float temp;
    float dist = 0;
    for (int j = 0; j < nrow; j++) {
        temp = norm(vec[j], origin);
        if (temp > dist) {
            dist = temp; start = j;
        }
    }

    id[0] = start; 
    int ref, temp_id; float ref_dis, temp_dis;
    vector<int> checked;
    checked.push_back(start);
    //cout << start << " " << vec[start][0] << " " << vec[start][1] << endl;
    for (int i = 1; i < nrow; i++) {
        ref = id[i-1]; 
        temp_dis = -1; temp_id;
        for (int j = 0; j < nrow; j++) {
            ref_dis = (ref < j)? mat[j][ref] : mat[ref][j];
            if ((ref_dis < temp_dis & ref != j) | temp_dis == -1) {
                if (notin(checked, j)) {
                    temp_dis = ref_dis; temp_id = j;
                }
            }
        }
        id[i] = temp_id; checked.push_back(temp_id);
        //cout << temp_id << " " << vec[temp_id][0] << " " << vec[temp_id][1] << endl;
    }
    
    for (int j = 0; j < nrow; j++) {
        id[j] = id[j] + 1;
    }
    return(id);
}