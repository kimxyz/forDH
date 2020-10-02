/******************************************************************************
 *                       Code generated with sympy 1.4                        *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_6613225906779492895);
void inv_err_fun(double *nom_x, double *true_x, double *out_116816818564585);
void H_mod_fun(double *state, double *out_5176169975441021617);
void f_fun(double *state, double dt, double *out_6693066138827208626);
void F_fun(double *state, double dt, double *out_5335318119803583076);
void h_3(double *state, double *unused, double *out_2297277348554745919);
void H_3(double *state, double *unused, double *out_3698690606063161161);
void h_4(double *state, double *unused, double *out_3418812439424126498);
void H_4(double *state, double *unused, double *out_8083655699941112067);
void h_9(double *state, double *unused, double *out_4587959295780742225);
void H_9(double *state, double *unused, double *out_331965737033253923);
void h_10(double *state, double *unused, double *out_2726123056400807542);
void H_10(double *state, double *unused, double *out_8367474854289666652);
void h_12(double *state, double *unused, double *out_6230486225429547787);
void H_12(double *state, double *unused, double *out_4435615913108327635);
void h_13(double *state, double *unused, double *out_3573558166983396831);
void H_13(double *state, double *unused, double *out_5168012595780274210);
void h_14(double *state, double *unused, double *out_4587959295780742225);
void H_14(double *state, double *unused, double *out_331965737033253923);
void h_19(double *state, double *unused, double *out_1972893041637311449);
void H_19(double *state, double *unused, double *out_5215323834937141177);
#define DIM 23
#define EDIM 22
#define MEDIM 22
typedef void (*Hfun)(double *, double *, double *);

void predict(double *x, double *P, double *Q, double dt);
const static double MAHA_THRESH_3 = 3.841459;
void update_3(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_4 = 7.814728;
void update_4(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_9 = 7.814728;
void update_9(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_10 = 7.814728;
void update_10(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_12 = 7.814728;
void update_12(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_13 = 7.814728;
void update_13(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_14 = 7.814728;
void update_14(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_19 = 7.814728;
void update_19(double *, double *, double *, double *, double *);