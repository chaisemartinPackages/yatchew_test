#include "stplugin.h"
#include <stdio.h>
#include <vector>
#include <cmath>
 
using namespace std;

void print(ST_double A) {
    char H[81];
    sprintf(H, "%f ", A);
    SF_display(H);
}

void nl() {
    char H[81];
    sprintf(H, "\n");
    SF_display(H);
}

bool notin(vector<ST_int> vec, ST_int num) {
    for (ST_int j = 0; j < vec.size(); j++) {
        if (vec[j] == num) {
            return(false);
        }
    }
    return(true);
}

ST_int min_e_dist(vector<vector<ST_double>> mat, ST_int p, vector<ST_int> vec) {
    ST_double dist = 0, temp; ST_int res = 0;
    for (ST_int j = 0; j < mat.size(); j++) {
        temp = (p - 1 > j)? mat[j][p-1] : mat[p-1][j];
        if (((dist != 0) & (temp < dist)) | (dist == 0)) {
            if (notin(vec, j+1)) {
                dist = temp; res = j+1;
            }
        }
    }
    return(res);
}

STDLL stata_call(int argc, char *argv[])
{   
    ST_int start, stop, nobs, nvars, S;
    ST_double value1, value2, d_temp;
    start = SF_in1();
    stop = SF_in2();
    nobs = stop - start + 1;
    nvars = SF_nvars()-1;

    vector<vector<ST_double>> M(nobs, vector<ST_double>(nobs));
    vector<ST_int> checked(nobs);

    for (ST_int i = start; i < start + nobs - 1; i++) {
        for (ST_int j = i + 1; j <= start + nobs - 1; j++) {
            d_temp = 0.0;
            for (ST_int c = 2; c <= nvars + 1; c++) {
                SF_vdata(c,i,&value1);SF_vdata(c,j,&value2);
                d_temp += pow((value1 - value2), 2.0);
            }
            M[i-start][j-start] = pow(d_temp, 0.5);
        }
    }

    SF_vstore(1, start, 1.0); S = 1.0;
    for (ST_int j = 2; j <= nobs; j++) {
        checked.push_back(S);
        S = min_e_dist(M, S, checked);
        SF_vstore(1, start + S - 1, j);
    }

    /*
    ST_int n_dis = nobs > 10? 10 : nobs;
    for (ST_int i = 0; i < n_dis; i ++) {
        for (ST_int j = 0; j < n_dis; j++ ) {
            print(M[i][j]);
        }
        nl();
    }
    */

    return((ST_retcode) 0);
}