/******************************************************************************
 *                       Code generated with sympy 1.4                        *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_2770769571783408795);
void inv_err_fun(double *nom_x, double *true_x, double *out_5284907962975002087);
void H_mod_fun(double *state, double *out_2132472091823993377);
void f_fun(double *state, double dt, double *out_8312258354532659232);
void F_fun(double *state, double dt, double *out_1195148538476340716);
void h_3(double *state, double *unused, double *out_122134289484665354);
void H_3(double *state, double *unused, double *out_267057401953282685);
void h_4(double *state, double *unused, double *out_8736897030297206565);
void H_4(double *state, double *unused, double *out_8312062215384217554);
void h_9(double *state, double *unused, double *out_5675586457956038187);
void H_9(double *state, double *unused, double *out_5726309947645159933);
void h_10(double *state, double *unused, double *out_6710352796085634921);
void H_10(double *state, double *unused, double *out_1195695374895787645);
void h_12(double *state, double *unused, double *out_2846239719562503159);
void H_12(double *state, double *unused, double *out_8361748188957713651);
void h_31(double *state, double *unused, double *out_5797297976849819690);
void H_31(double *state, double *unused, double *out_8455250968147141189);
void h_32(double *state, double *unused, double *out_6815683808841322701);
void H_32(double *state, double *unused, double *out_1531751758097609597);
void h_13(double *state, double *unused, double *out_8900641810197606043);
void H_13(double *state, double *unused, double *out_8039295433079851980);
void h_14(double *state, double *unused, double *out_5675586457956038187);
void H_14(double *state, double *unused, double *out_5726309947645159933);
void h_19(double *state, double *unused, double *out_1057366158597397627);
void H_19(double *state, double *unused, double *out_9088250110993165863);
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
const static double MAHA_THRESH_31 = 7.814728;
void update_31(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_32 = 9.487729;
void update_32(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_13 = 7.814728;
void update_13(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_14 = 7.814728;
void update_14(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_19 = 7.814728;
void update_19(double *, double *, double *, double *, double *);