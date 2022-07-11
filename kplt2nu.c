#include <stdio.h>
#include <math.h>
#include "bisnwt.h"

#define PI 3.1415926535897932384626

struct kepler_prm {
    double e, M;
};

double amitj(double t, double t_p, double T);
double eq_kepler(double E, void* prm);
double fkepler(double E, void* prm);
double get_v(double E, double e);


int main (int argc, char *argv[]) {
    double e, T, M0, tf;
    int nt;

    if (argc < 6
        || sscanf(argv[1], "%lf", &e) != 1
        || sscanf(argv[2], "%lf", &T) != 1
        || sscanf(argv[3], "%lf", &M0) != 1
        || sscanf(argv[4], "%lf", &tf) != 1
        || sscanf(argv[5], "%d", &nt) != 1
            ) {
        fprintf(stderr, "%s e T M0 tf nt\n", argv[0]);
        return -1;
    }

    double dlt = 2.5, maxit = 10, tol = 1e-12;
    double t = 0, tp = - (T * M0) / (2 * PI);
    double a, b, E = 0;
    double M, v;
    int k;
    double vect[2];

    vect[0] = e; //en aquest vector guardarem els valors e i posteriorment l'anomalia mitjana

    for (int i = 0; i <= nt; i++) {
        M = amitj(t, tp, T);
        vect[1] = M;
        void *prm = vect;
        a = M - PI;  //inici interval
        b = M + PI * PI;  //final interval
        bisnwt(a, b, &E, &dlt, tol, maxit, &eq_kepler, &fkepler, prm);
        k = floor(E / (2 * PI));
        v = get_v(E, e);
        //triem la determinació de l'arcos
        if (E <= k * 2 * PI + PI) {
            v += k * 2 * PI;
        } else {
            v = 2 * PI - v + k * 2 * PI;
        }

        printf("%.16G ",t);
        printf("%.16G ",M);
        printf("%.16G \n",v);

        t += tf / nt;
    }
    return 0;
}



double amitj(double t, double tp, double T){
    return ((t - tp) / T) * 2 * PI;
}

double eq_kepler (double E, void *prm) {
    struct kepler_prm *pkprm=prm;
    return E-(pkprm->e)*sin(E)-(pkprm->M);
}

double fkepler(double E, void* prm){
    struct kepler_prm *pkprm=prm;
    return 1 - pkprm->e * cos(E);
}

double get_v(double E, double e){
    return acos((e - cos(E)) / (e * cos(E) - 1));
}