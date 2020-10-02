/******************************************************************************
 *                       Code generated with sympy 1.4                        *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_2970325197894364982);
void inv_err_fun(double *nom_x, double *true_x, double *out_2772786700340426778);
void H_mod_fun(double *state, double *out_3535125193612146439);
void f_fun(double *state, double dt, double *out_9101139339523934526);
void F_fun(double *state, double dt, double *out_1730133153759752138);
void h_25(double *state, double *unused, double *out_2281010244982640375);
void H_25(double *state, double *unused, double *out_8190463851841531169);
void h_24(double *state, double *unused, double *out_6158622900840544827);
void H_24(double *state, double *unused, double *out_1318428548392123251);
void h_30(double *state, double *unused, double *out_6699036697577757067);
void H_30(double *state, double *unused, double *out_3312134217808160971);
void h_26(double *state, double *unused, double *out_464548907851401214);
void H_26(double *state, double *unused, double *out_8990925556697705385);
void h_27(double *state, double *unused, double *out_6149024344084742650);
void H_27(double *state, double *unused, double *out_6776346510559708411);
void h_29(double *state, double *unused, double *out_3948524039439182520);
void H_29(double *state, double *unused, double *out_73535388410856529);
void h_28(double *state, double *unused, double *out_57999947746821025);
void H_28(double *state, double *unused, double *out_6153668077513873847);
#define DIM 8
#define EDIM 8
#define MEDIM 8
typedef void (*Hfun)(double *, double *, double *);

void predict(double *x, double *P, double *Q, double dt);
const static double MAHA_THRESH_25 = 3.841459;
void update_25(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_24 = 5.991465;
void update_24(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_30 = 3.841459;
void update_30(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_26 = 3.841459;
void update_26(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_27 = 3.841459;
void update_27(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_29 = 3.841459;
void update_29(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_28 = 5.991465;
void update_28(double *, double *, double *, double *, double *);
void set_mass(double x);

void set_rotational_inertia(double x);

void set_center_to_front(double x);

void set_center_to_rear(double x);

void set_stiffness_front(double x);

void set_stiffness_rear(double x);
