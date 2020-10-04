/******************************************************************************
 *                       Code generated with sympy 1.4                        *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_8597576853528407995);
void inv_err_fun(double *nom_x, double *true_x, double *out_337191706622251325);
void H_mod_fun(double *state, double *out_7252594450031301036);
void f_fun(double *state, double dt, double *out_750442287767465423);
void F_fun(double *state, double dt, double *out_8768775507318021386);
void h_3(double *state, double *unused, double *out_1402387130017346519);
void H_3(double *state, double *unused, double *out_6155385720429878687);
void h_4(double *state, double *unused, double *out_7748716648231134878);
void H_4(double *state, double *unused, double *out_1517436953518315750);
void h_9(double *state, double *unused, double *out_7532412872684409159);
void H_9(double *state, double *unused, double *out_2432857147054178634);
void h_10(double *state, double *unused, double *out_6275582451076850556);
void H_10(double *state, double *unused, double *out_1671852089301200390);
void h_12(double *state, double *unused, double *out_8201598838486474843);
void H_12(double *state, double *unused, double *out_6120524722267320742);
void h_13(double *state, double *unused, double *out_331065800735883660);
void H_13(double *state, double *unused, double *out_6494004337227518786);
void h_14(double *state, double *unused, double *out_7532412872684409159);
void H_14(double *state, double *unused, double *out_2432857147054178634);
void h_19(double *state, double *unused, double *out_7293303393803558751);
void H_19(double *state, double *unused, double *out_2543398202372958174);
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