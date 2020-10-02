
extern "C"{

double mass;

void set_mass(double x){ mass = x;}

double rotational_inertia;

void set_rotational_inertia(double x){ rotational_inertia = x;}

double center_to_front;

void set_center_to_front(double x){ center_to_front = x;}

double center_to_rear;

void set_center_to_rear(double x){ center_to_rear = x;}

double stiffness_front;

void set_stiffness_front(double x){ stiffness_front = x;}

double stiffness_rear;

void set_stiffness_rear(double x){ stiffness_rear = x;}

}
extern "C" {
#include <math.h>
/******************************************************************************
 *                       Code generated with sympy 1.4                        *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_2970325197894364982) {
   out_2970325197894364982[0] = delta_x[0] + nom_x[0];
   out_2970325197894364982[1] = delta_x[1] + nom_x[1];
   out_2970325197894364982[2] = delta_x[2] + nom_x[2];
   out_2970325197894364982[3] = delta_x[3] + nom_x[3];
   out_2970325197894364982[4] = delta_x[4] + nom_x[4];
   out_2970325197894364982[5] = delta_x[5] + nom_x[5];
   out_2970325197894364982[6] = delta_x[6] + nom_x[6];
   out_2970325197894364982[7] = delta_x[7] + nom_x[7];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_2772786700340426778) {
   out_2772786700340426778[0] = -nom_x[0] + true_x[0];
   out_2772786700340426778[1] = -nom_x[1] + true_x[1];
   out_2772786700340426778[2] = -nom_x[2] + true_x[2];
   out_2772786700340426778[3] = -nom_x[3] + true_x[3];
   out_2772786700340426778[4] = -nom_x[4] + true_x[4];
   out_2772786700340426778[5] = -nom_x[5] + true_x[5];
   out_2772786700340426778[6] = -nom_x[6] + true_x[6];
   out_2772786700340426778[7] = -nom_x[7] + true_x[7];
}
void H_mod_fun(double *state, double *out_3535125193612146439) {
   out_3535125193612146439[0] = 1.0;
   out_3535125193612146439[1] = 0.0;
   out_3535125193612146439[2] = 0.0;
   out_3535125193612146439[3] = 0.0;
   out_3535125193612146439[4] = 0.0;
   out_3535125193612146439[5] = 0.0;
   out_3535125193612146439[6] = 0.0;
   out_3535125193612146439[7] = 0.0;
   out_3535125193612146439[8] = 0.0;
   out_3535125193612146439[9] = 1.0;
   out_3535125193612146439[10] = 0.0;
   out_3535125193612146439[11] = 0.0;
   out_3535125193612146439[12] = 0.0;
   out_3535125193612146439[13] = 0.0;
   out_3535125193612146439[14] = 0.0;
   out_3535125193612146439[15] = 0.0;
   out_3535125193612146439[16] = 0.0;
   out_3535125193612146439[17] = 0.0;
   out_3535125193612146439[18] = 1.0;
   out_3535125193612146439[19] = 0.0;
   out_3535125193612146439[20] = 0.0;
   out_3535125193612146439[21] = 0.0;
   out_3535125193612146439[22] = 0.0;
   out_3535125193612146439[23] = 0.0;
   out_3535125193612146439[24] = 0.0;
   out_3535125193612146439[25] = 0.0;
   out_3535125193612146439[26] = 0.0;
   out_3535125193612146439[27] = 1.0;
   out_3535125193612146439[28] = 0.0;
   out_3535125193612146439[29] = 0.0;
   out_3535125193612146439[30] = 0.0;
   out_3535125193612146439[31] = 0.0;
   out_3535125193612146439[32] = 0.0;
   out_3535125193612146439[33] = 0.0;
   out_3535125193612146439[34] = 0.0;
   out_3535125193612146439[35] = 0.0;
   out_3535125193612146439[36] = 1.0;
   out_3535125193612146439[37] = 0.0;
   out_3535125193612146439[38] = 0.0;
   out_3535125193612146439[39] = 0.0;
   out_3535125193612146439[40] = 0.0;
   out_3535125193612146439[41] = 0.0;
   out_3535125193612146439[42] = 0.0;
   out_3535125193612146439[43] = 0.0;
   out_3535125193612146439[44] = 0.0;
   out_3535125193612146439[45] = 1.0;
   out_3535125193612146439[46] = 0.0;
   out_3535125193612146439[47] = 0.0;
   out_3535125193612146439[48] = 0.0;
   out_3535125193612146439[49] = 0.0;
   out_3535125193612146439[50] = 0.0;
   out_3535125193612146439[51] = 0.0;
   out_3535125193612146439[52] = 0.0;
   out_3535125193612146439[53] = 0.0;
   out_3535125193612146439[54] = 1.0;
   out_3535125193612146439[55] = 0.0;
   out_3535125193612146439[56] = 0.0;
   out_3535125193612146439[57] = 0.0;
   out_3535125193612146439[58] = 0.0;
   out_3535125193612146439[59] = 0.0;
   out_3535125193612146439[60] = 0.0;
   out_3535125193612146439[61] = 0.0;
   out_3535125193612146439[62] = 0.0;
   out_3535125193612146439[63] = 1.0;
}
void f_fun(double *state, double dt, double *out_9101139339523934526) {
   out_9101139339523934526[0] = state[0];
   out_9101139339523934526[1] = state[1];
   out_9101139339523934526[2] = state[2];
   out_9101139339523934526[3] = state[3];
   out_9101139339523934526[4] = state[4];
   out_9101139339523934526[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_9101139339523934526[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_9101139339523934526[7] = state[7];
}
void F_fun(double *state, double dt, double *out_1730133153759752138) {
   out_1730133153759752138[0] = 1;
   out_1730133153759752138[1] = 0;
   out_1730133153759752138[2] = 0;
   out_1730133153759752138[3] = 0;
   out_1730133153759752138[4] = 0;
   out_1730133153759752138[5] = 0;
   out_1730133153759752138[6] = 0;
   out_1730133153759752138[7] = 0;
   out_1730133153759752138[8] = 0;
   out_1730133153759752138[9] = 1;
   out_1730133153759752138[10] = 0;
   out_1730133153759752138[11] = 0;
   out_1730133153759752138[12] = 0;
   out_1730133153759752138[13] = 0;
   out_1730133153759752138[14] = 0;
   out_1730133153759752138[15] = 0;
   out_1730133153759752138[16] = 0;
   out_1730133153759752138[17] = 0;
   out_1730133153759752138[18] = 1;
   out_1730133153759752138[19] = 0;
   out_1730133153759752138[20] = 0;
   out_1730133153759752138[21] = 0;
   out_1730133153759752138[22] = 0;
   out_1730133153759752138[23] = 0;
   out_1730133153759752138[24] = 0;
   out_1730133153759752138[25] = 0;
   out_1730133153759752138[26] = 0;
   out_1730133153759752138[27] = 1;
   out_1730133153759752138[28] = 0;
   out_1730133153759752138[29] = 0;
   out_1730133153759752138[30] = 0;
   out_1730133153759752138[31] = 0;
   out_1730133153759752138[32] = 0;
   out_1730133153759752138[33] = 0;
   out_1730133153759752138[34] = 0;
   out_1730133153759752138[35] = 0;
   out_1730133153759752138[36] = 1;
   out_1730133153759752138[37] = 0;
   out_1730133153759752138[38] = 0;
   out_1730133153759752138[39] = 0;
   out_1730133153759752138[40] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_1730133153759752138[41] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_1730133153759752138[42] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_1730133153759752138[43] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_1730133153759752138[44] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_1730133153759752138[45] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_1730133153759752138[46] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_1730133153759752138[47] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_1730133153759752138[48] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_1730133153759752138[49] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_1730133153759752138[50] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_1730133153759752138[51] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_1730133153759752138[52] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_1730133153759752138[53] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_1730133153759752138[54] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_1730133153759752138[55] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_1730133153759752138[56] = 0;
   out_1730133153759752138[57] = 0;
   out_1730133153759752138[58] = 0;
   out_1730133153759752138[59] = 0;
   out_1730133153759752138[60] = 0;
   out_1730133153759752138[61] = 0;
   out_1730133153759752138[62] = 0;
   out_1730133153759752138[63] = 1;
}
void h_25(double *state, double *unused, double *out_2281010244982640375) {
   out_2281010244982640375[0] = state[6];
}
void H_25(double *state, double *unused, double *out_8190463851841531169) {
   out_8190463851841531169[0] = 0;
   out_8190463851841531169[1] = 0;
   out_8190463851841531169[2] = 0;
   out_8190463851841531169[3] = 0;
   out_8190463851841531169[4] = 0;
   out_8190463851841531169[5] = 0;
   out_8190463851841531169[6] = 1;
   out_8190463851841531169[7] = 0;
}
void h_24(double *state, double *unused, double *out_6158622900840544827) {
   out_6158622900840544827[0] = state[4];
   out_6158622900840544827[1] = state[5];
}
void H_24(double *state, double *unused, double *out_1318428548392123251) {
   out_1318428548392123251[0] = 0;
   out_1318428548392123251[1] = 0;
   out_1318428548392123251[2] = 0;
   out_1318428548392123251[3] = 0;
   out_1318428548392123251[4] = 1;
   out_1318428548392123251[5] = 0;
   out_1318428548392123251[6] = 0;
   out_1318428548392123251[7] = 0;
   out_1318428548392123251[8] = 0;
   out_1318428548392123251[9] = 0;
   out_1318428548392123251[10] = 0;
   out_1318428548392123251[11] = 0;
   out_1318428548392123251[12] = 0;
   out_1318428548392123251[13] = 1;
   out_1318428548392123251[14] = 0;
   out_1318428548392123251[15] = 0;
}
void h_30(double *state, double *unused, double *out_6699036697577757067) {
   out_6699036697577757067[0] = state[4];
}
void H_30(double *state, double *unused, double *out_3312134217808160971) {
   out_3312134217808160971[0] = 0;
   out_3312134217808160971[1] = 0;
   out_3312134217808160971[2] = 0;
   out_3312134217808160971[3] = 0;
   out_3312134217808160971[4] = 1;
   out_3312134217808160971[5] = 0;
   out_3312134217808160971[6] = 0;
   out_3312134217808160971[7] = 0;
}
void h_26(double *state, double *unused, double *out_464548907851401214) {
   out_464548907851401214[0] = state[7];
}
void H_26(double *state, double *unused, double *out_8990925556697705385) {
   out_8990925556697705385[0] = 0;
   out_8990925556697705385[1] = 0;
   out_8990925556697705385[2] = 0;
   out_8990925556697705385[3] = 0;
   out_8990925556697705385[4] = 0;
   out_8990925556697705385[5] = 0;
   out_8990925556697705385[6] = 0;
   out_8990925556697705385[7] = 1;
}
void h_27(double *state, double *unused, double *out_6149024344084742650) {
   out_6149024344084742650[0] = state[3];
}
void H_27(double *state, double *unused, double *out_6776346510559708411) {
   out_6776346510559708411[0] = 0;
   out_6776346510559708411[1] = 0;
   out_6776346510559708411[2] = 0;
   out_6776346510559708411[3] = 1;
   out_6776346510559708411[4] = 0;
   out_6776346510559708411[5] = 0;
   out_6776346510559708411[6] = 0;
   out_6776346510559708411[7] = 0;
}
void h_29(double *state, double *unused, double *out_3948524039439182520) {
   out_3948524039439182520[0] = state[1];
}
void H_29(double *state, double *unused, double *out_73535388410856529) {
   out_73535388410856529[0] = 0;
   out_73535388410856529[1] = 1;
   out_73535388410856529[2] = 0;
   out_73535388410856529[3] = 0;
   out_73535388410856529[4] = 0;
   out_73535388410856529[5] = 0;
   out_73535388410856529[6] = 0;
   out_73535388410856529[7] = 0;
}
void h_28(double *state, double *unused, double *out_57999947746821025) {
   out_57999947746821025[0] = state[5];
   out_57999947746821025[1] = state[6];
}
void H_28(double *state, double *unused, double *out_6153668077513873847) {
   out_6153668077513873847[0] = 0;
   out_6153668077513873847[1] = 0;
   out_6153668077513873847[2] = 0;
   out_6153668077513873847[3] = 0;
   out_6153668077513873847[4] = 0;
   out_6153668077513873847[5] = 1;
   out_6153668077513873847[6] = 0;
   out_6153668077513873847[7] = 0;
   out_6153668077513873847[8] = 0;
   out_6153668077513873847[9] = 0;
   out_6153668077513873847[10] = 0;
   out_6153668077513873847[11] = 0;
   out_6153668077513873847[12] = 0;
   out_6153668077513873847[13] = 0;
   out_6153668077513873847[14] = 1;
   out_6153668077513873847[15] = 0;
}
}

extern "C"{
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
}

#include <eigen3/Eigen/Dense>
#include <iostream>

typedef Eigen::Matrix<double, DIM, DIM, Eigen::RowMajor> DDM;
typedef Eigen::Matrix<double, EDIM, EDIM, Eigen::RowMajor> EEM;
typedef Eigen::Matrix<double, DIM, EDIM, Eigen::RowMajor> DEM;

void predict(double *in_x, double *in_P, double *in_Q, double dt) {
  typedef Eigen::Matrix<double, MEDIM, MEDIM, Eigen::RowMajor> RRM;
  
  double nx[DIM] = {0};
  double in_F[EDIM*EDIM] = {0};

  // functions from sympy
  f_fun(in_x, dt, nx);
  F_fun(in_x, dt, in_F);


  EEM F(in_F);
  EEM P(in_P);
  EEM Q(in_Q);

  RRM F_main = F.topLeftCorner(MEDIM, MEDIM);
  P.topLeftCorner(MEDIM, MEDIM) = (F_main * P.topLeftCorner(MEDIM, MEDIM)) * F_main.transpose();
  P.topRightCorner(MEDIM, EDIM - MEDIM) = F_main * P.topRightCorner(MEDIM, EDIM - MEDIM);
  P.bottomLeftCorner(EDIM - MEDIM, MEDIM) = P.bottomLeftCorner(EDIM - MEDIM, MEDIM) * F_main.transpose();

  P = P + dt*Q;

  // copy out state
  memcpy(in_x, nx, DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
}

// note: extra_args dim only correct when null space projecting
// otherwise 1
template <int ZDIM, int EADIM, bool MAHA_TEST>
void update(double *in_x, double *in_P, Hfun h_fun, Hfun H_fun, Hfun Hea_fun, double *in_z, double *in_R, double *in_ea, double MAHA_THRESHOLD) {
  typedef Eigen::Matrix<double, ZDIM, ZDIM, Eigen::RowMajor> ZZM;
  typedef Eigen::Matrix<double, ZDIM, DIM, Eigen::RowMajor> ZDM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, EDIM, Eigen::RowMajor> XEM;
  //typedef Eigen::Matrix<double, EDIM, ZDIM, Eigen::RowMajor> EZM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, 1> X1M;
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> XXM;

  double in_hx[ZDIM] = {0};
  double in_H[ZDIM * DIM] = {0};
  double in_H_mod[EDIM * DIM] = {0};
  double delta_x[EDIM] = {0};
  double x_new[DIM] = {0};


  // state x, P
  Eigen::Matrix<double, ZDIM, 1> z(in_z);
  EEM P(in_P);
  ZZM pre_R(in_R);
  
  // functions from sympy
  h_fun(in_x, in_ea, in_hx);
  H_fun(in_x, in_ea, in_H);
  ZDM pre_H(in_H); 
  
  // get y (y = z - hx)
  Eigen::Matrix<double, ZDIM, 1> pre_y(in_hx); pre_y = z - pre_y;
  X1M y; XXM H; XXM R;
  if (Hea_fun){
    typedef Eigen::Matrix<double, ZDIM, EADIM, Eigen::RowMajor> ZAM;
    double in_Hea[ZDIM * EADIM] = {0};
    Hea_fun(in_x, in_ea, in_Hea);
    ZAM Hea(in_Hea);
    XXM A = Hea.transpose().fullPivLu().kernel();
   

    y = A.transpose() * pre_y;
    H = A.transpose() * pre_H;
    R = A.transpose() * pre_R * A;
  } else {
    y = pre_y;
    H = pre_H;
    R = pre_R;
  }
  // get modified H
  H_mod_fun(in_x, in_H_mod);
  DEM H_mod(in_H_mod);
  XEM H_err = H * H_mod;
  
  // Do mahalobis distance test
  if (MAHA_TEST){
    XXM a = (H_err * P * H_err.transpose() + R).inverse();
    double maha_dist = y.transpose() * a * y;
    if (maha_dist > MAHA_THRESHOLD){
      R = 1.0e16 * R;
    }
  }

  // Outlier resilient weighting
  double weight = 1;//(1.5)/(1 + y.squaredNorm()/R.sum());

  // kalman gains and I_KH
  XXM S = ((H_err * P) * H_err.transpose()) + R/weight;
  XEM KT = S.fullPivLu().solve(H_err * P.transpose());
  //EZM K = KT.transpose(); TODO: WHY DOES THIS NOT COMPILE?
  //EZM K = S.fullPivLu().solve(H_err * P.transpose()).transpose();
  //std::cout << "Here is the matrix rot:\n" << K << std::endl;
  EEM I_KH = Eigen::Matrix<double, EDIM, EDIM>::Identity() - (KT.transpose() * H_err);

  // update state by injecting dx
  Eigen::Matrix<double, EDIM, 1> dx(delta_x);
  dx  = (KT.transpose() * y);
  memcpy(delta_x, dx.data(), EDIM * sizeof(double));
  err_fun(in_x, delta_x, x_new);
  Eigen::Matrix<double, DIM, 1> x(x_new);
 
  // update cov 
  P = ((I_KH * P) * I_KH.transpose()) + ((KT.transpose() * R) * KT);

  // copy out state
  memcpy(in_x, x.data(), DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
  memcpy(in_z, y.data(), y.rows() * sizeof(double));
}



extern "C"{

      void update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_25, H_25, NULL, in_z, in_R, in_ea, MAHA_THRESH_25);
      }
    
      void update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<2,3,0>(in_x, in_P, h_24, H_24, NULL, in_z, in_R, in_ea, MAHA_THRESH_24);
      }
    
      void update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_30, H_30, NULL, in_z, in_R, in_ea, MAHA_THRESH_30);
      }
    
      void update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_26, H_26, NULL, in_z, in_R, in_ea, MAHA_THRESH_26);
      }
    
      void update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_27, H_27, NULL, in_z, in_R, in_ea, MAHA_THRESH_27);
      }
    
      void update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_29, H_29, NULL, in_z, in_R, in_ea, MAHA_THRESH_29);
      }
    
      void update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<2,3,0>(in_x, in_P, h_28, H_28, NULL, in_z, in_R, in_ea, MAHA_THRESH_28);
      }
    
}
