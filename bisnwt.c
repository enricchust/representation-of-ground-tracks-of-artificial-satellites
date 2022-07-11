#include <stdio.h>
#include <math.h>
#include <string.h>
#include "bisnwt.h"


int bisnwt (double a, double b, double *arr,
            double *dlt, double tol, int maxit,
            double (*f)(double,void*), double (*df)(double,void*), void *prm){
    double fa, fc;
    double c = (a + b) / 2; //punt mig
    int iteracions = 0;

    fa = (*f)(a,prm);

    while(b - a > *dlt){
        fc = (*f)(c, prm);
        if(fa * fc > 0){  //mirem si hi ha o no un zero
            a = c;
            fa = fc;
        } else if (fa * fc < 0){
            b = c;
        } else break;
        c = (a + b) / 2;
    }

    if(*dlt <= tol){
        *arr = c;
        return -1;
    }

    //apliquem el mètode Newton
    double xn = c, xn1 = 1;
    while(fabs(xn - xn1) >= tol && iteracions < maxit){
        xn1 = xn;
        xn -= (*f)(xn, prm) / (*df)(xn, prm);
        iteracions++;
    }
    //si no convergeix
    if(fabs(xn - xn1) >= tol){
        *dlt = *dlt / 2;
        return bisnwt(a, b, arr, dlt, tol, maxit,  f,  df, prm);
    }else{
        //si convergeix
        *arr = xn;
        return iteracions;
    }
}