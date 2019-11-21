#include "VisionSolverMPC.h"
#include "../convexMPC/common_types.h"
#include "VisionMPC_interface.h"
#include "VisionRobotState.h"
#include <eigen3/Eigen/Dense>
#include <cmath>
#include <eigen3/unsupported/Eigen/MatrixFunctions>
#include <qpOASES.hpp>
#include <stdio.h>
#include <sys/time.h>

#define V_BIG_NUMBER 5e10
//big enough to act like infinity, small enough to avoid numerical weirdness.

VisionRobotState v_rs;
using std::cout;
using std::endl;
using Eigen::Dynamic;

Matrix<fpt,Dynamic,13> vA_qp;
Matrix<fpt,Dynamic,Dynamic> vB_qp;
Matrix<fpt,13,12> vBdt;
Matrix<fpt,13,13> vAdt;
Matrix<fpt,25,25> vABc,v_expmm;
Matrix<fpt,Dynamic,Dynamic> vS;
Matrix<fpt,Dynamic,1> vX_d;
Matrix<fpt,Dynamic,1> vU_b;
Matrix<fpt,Dynamic,Dynamic> v_fmat;

Matrix<fpt,Dynamic,Dynamic> v_qH;
Matrix<fpt,Dynamic,1> v_qg;

Matrix<fpt,Dynamic,Dynamic> v_eye_12h;

qpOASES::real_t* vH_qpoases;
qpOASES::real_t* vg_qpoases;
qpOASES::real_t* vA_qpoases;
qpOASES::real_t* vlb_qpoases;
qpOASES::real_t* vub_qpoases;
qpOASES::real_t* vq_soln;

qpOASES::real_t* vH_red;
qpOASES::real_t* vg_red;
qpOASES::real_t* vA_red;
qpOASES::real_t* vlb_red;
qpOASES::real_t* vub_red;
qpOASES::real_t* vq_red;
u8 v_real_allocated = 0;


char v_var_elim[2000];
char v_con_elim[2000];

mfp* vision_get_q_soln() {
  return vq_soln;
}

s8 v_near_zero(fpt a)
{
  return (a < 0.01 && a > -.01) ;
}

s8 v_near_one(fpt a)
{
  return v_near_zero(a-1);
}
void v_matrix_to_real(qpOASES::real_t* dst, Matrix<fpt,Dynamic,Dynamic> src, s16 rows, s16 cols)
{
  s32 a = 0;
  for(s16 r = 0; r < rows; r++)
  {
    for(s16 c = 0; c < cols; c++)
    {
      dst[a] = src(r,c);
      a++;
    }
  }
}


void vision_c2qp(Matrix<fpt,13,13> Ac, Matrix<fpt,13,12> Bc,fpt dt,s16 horizon)
{
  vABc.setZero();
  vABc.block(0,0,13,13) = Ac;
  vABc.block(0,13,13,12) = Bc;
  vABc = dt*vABc;
  v_expmm = vABc.exp();
  vAdt = v_expmm.block(0,0,13,13);
  vBdt = v_expmm.block(0,13,13,12);
  if(horizon > 19) {
    throw std::runtime_error("horizon is too long!");
  }

  Matrix<fpt,13,13> powerMats[20];
  powerMats[0].setIdentity();
  for(int i = 1; i < horizon+1; i++) {
    powerMats[i] = vAdt * powerMats[i-1];
  }

  for(s16 r = 0; r < horizon; r++)
  {
    vA_qp.block(13*r,0,13,13) = powerMats[r+1];
    for(s16 c = 0; c < horizon; c++)
    {
      if(r >= c)
      {
        s16 a_num = r-c;
        vB_qp.block(13*r,12*c,13,12) = powerMats[a_num] * vBdt;
      }
    }
  }

}

void vision_resize_qp_mats(s16 horizon)
{
  int mcount = 0;
  int h2 = horizon*horizon;

  vA_qp.resize(13*horizon, Eigen::NoChange);
  mcount += 13*horizon*1;

  vB_qp.resize(13*horizon, 12*horizon);
  mcount += 13*h2*12;

  vS.resize(13*horizon, 13*horizon);
  mcount += 13*13*h2;

  vX_d.resize(13*horizon, Eigen::NoChange);
  mcount += 13*horizon;

  vU_b.resize(20*horizon, Eigen::NoChange);
  mcount += 20*horizon;

  v_fmat.resize(20*horizon, 12*horizon);
  mcount += 20*12*h2;

  v_qH.resize(12*horizon, 12*horizon);
  mcount += 12*12*h2;

  v_qg.resize(12*horizon, Eigen::NoChange);
  mcount += 12*horizon;

  v_eye_12h.resize(12*horizon, 12*horizon);
  mcount += 12*12*horizon;

  //printf("realloc'd %d floating point numbers.\n",mcount);
  mcount = 0;

  vA_qp.setZero();
  vB_qp.setZero();
  vS.setZero();
  vX_d.setZero();
  vU_b.setZero();
  v_fmat.setZero();
  v_qH.setZero();
  v_eye_12h.setIdentity();

  //TODO: use realloc instead of free/malloc on size changes

  if(v_real_allocated)
  {

    free(vH_qpoases);
    free(vg_qpoases);
    free(vA_qpoases);
    free(vlb_qpoases);
    free(vub_qpoases);
    free(vq_soln);
    free(vH_red);
    free(vg_red);
    free(vA_red);
    free(vlb_red);
    free(vub_red);
    free(vq_red);
  }

  vH_qpoases = (qpOASES::real_t*)malloc(12*12*horizon*horizon*sizeof(qpOASES::real_t));
  mcount += 12*12*h2;
  vg_qpoases = (qpOASES::real_t*)malloc(12*1*horizon*sizeof(qpOASES::real_t));
  mcount += 12*horizon;
  vA_qpoases = (qpOASES::real_t*)malloc(12*20*horizon*horizon*sizeof(qpOASES::real_t));
  mcount += 12*20*h2;
  vlb_qpoases = (qpOASES::real_t*)malloc(20*1*horizon*sizeof(qpOASES::real_t));
  mcount += 20*horizon;
  vub_qpoases = (qpOASES::real_t*)malloc(20*1*horizon*sizeof(qpOASES::real_t));
  mcount += 20*horizon;
  vq_soln = (qpOASES::real_t*)malloc(12*horizon*sizeof(qpOASES::real_t));
  mcount += 12*horizon;

  vH_red = (qpOASES::real_t*)malloc(12*12*horizon*horizon*sizeof(qpOASES::real_t));
  mcount += 12*12*h2;
  vg_red = (qpOASES::real_t*)malloc(12*1*horizon*sizeof(qpOASES::real_t));
  mcount += 12*horizon;
  vA_red = (qpOASES::real_t*)malloc(12*20*horizon*horizon*sizeof(qpOASES::real_t));
  mcount += 12*20*h2;
  vlb_red = (qpOASES::real_t*)malloc(20*1*horizon*sizeof(qpOASES::real_t));
  mcount += 20*horizon;
  vub_red = (qpOASES::real_t*)malloc(20*1*horizon*sizeof(qpOASES::real_t));
  mcount += 20*horizon;
  vq_red = (qpOASES::real_t*)malloc(12*horizon*sizeof(qpOASES::real_t));
  mcount += 12*horizon;
  v_real_allocated = 1;

  //printf("malloc'd %d floating point numbers.\n",mcount);



#ifdef K_DEBUG
  printf("RESIZED MATRICES FOR HORIZON: %d\n",horizon);
#endif
}

inline Matrix<fpt,3,3> cross_mat(Matrix<fpt,3,3> I_inv, Matrix<fpt,3,1> r)
{
  Matrix<fpt,3,3> cm;
  cm << 0.f, -r(2), r(1),
     r(2), 0.f, -r(0),
     -r(1), r(0), 0.f;
  return I_inv * cm;
}
//continuous time state space matrices.
void vision_ct_ss_mats(Matrix<fpt,3,3> vI_world, fpt m, Matrix<fpt,3,4> r_feet, 
    Matrix<fpt,3,3> R_yaw, Matrix<fpt,13,13>& A, Matrix<fpt,13,12>& B, float x_drag)
{
  A.setZero();
  A(3,9) = 1.f;
  A(9,9) = x_drag;
  A(4,10) = 1.f;
  A(5,11) = 1.f;

  A(11,12) = 1.f;
  A.block(0,6,3,3) = R_yaw.transpose();

  B.setZero();
  Matrix<fpt,3,3> I_inv = vI_world.inverse();

  for(s16 b = 0; b < 4; b++)
  {
    B.block(6,b*3,3,3) = cross_mat(I_inv,r_feet.col(b));
    B.block(9,b*3,3,3) = Matrix<fpt,3,3>::Identity() / m;
  }
}


void vision_quat_to_rpy(Quaternionf q, Matrix<fpt,3,1>& rpy)
{
  //from my MATLAB implementation

  //edge case!
  fpt as = vision_t_min(-2.*(q.x()*q.z()-q.w()*q.y()),.99999);
  rpy(0) = atan2(2.f*(q.x()*q.y()+q.w()*q.z()),vision_sq(q.w()) + vision_sq(q.x()) - vision_sq(q.y()) - vision_sq(q.z()));
  rpy(1) = asin(as);
  rpy(2) = atan2(2.f*(q.y()*q.z()+q.w()*q.x()),vision_sq(q.w()) - vision_sq(q.x()) - vision_sq(q.y()) + vision_sq(q.z()));

}

Matrix<fpt,13,1> v_x_0;
Matrix<fpt,3,3> vI_world;
Matrix<fpt,13,13> vA_ct;
Matrix<fpt,13,12> vB_ct_r;


void vision_solve_mpc(vision_mpc_update_data_t* update, vision_mpc_problem_setup* setup)
{
  v_rs.set(update->p, update->v, update->q, update->w, update->r, update->yaw);

  //roll pitch yaw
  Matrix<fpt,3,1> rpy;
  vision_quat_to_rpy(v_rs.q,rpy);

  //initial state (13 state representation)
  v_x_0 << rpy(2), rpy(1), rpy(0), v_rs.p , v_rs.w, v_rs.v, -9.8f;
  vI_world = v_rs.R_yaw * v_rs.I_body * v_rs.R_yaw.transpose(); //original
  vision_ct_ss_mats(vI_world,v_rs.m,v_rs.r_feet,v_rs.R_yaw,vA_ct,vB_ct_r, update->x_drag);


  //QP matrices
  vision_c2qp(vA_ct,vB_ct_r,setup->dt,setup->horizon);

  //weights
  Matrix<fpt,13,1> full_weight;
  for(u8 i = 0; i < 12; i++)
    full_weight(i) = update->weights[i];
  full_weight(12) = 0.f;
  vS.diagonal() = full_weight.replicate(setup->horizon,1);

  //trajectory
  for(s16 i = 0; i < setup->horizon; i++)
  {
    for(s16 j = 0; j < 12; j++)
      vX_d(13*i+j,0) = update->traj[12*i+j];
  }
  //cout<<"XD:\n"<<vX_d<<endl;



  //note - I'm not doing the shifting here.
  s16 k = 0;
  for(s16 i = 0; i < setup->horizon; i++)
  {
    for(s16 j = 0; j < 4; j++)
    {
      vU_b(5*k + 0) = V_BIG_NUMBER;
      vU_b(5*k + 1) = V_BIG_NUMBER;
      vU_b(5*k + 2) = V_BIG_NUMBER;
      vU_b(5*k + 3) = V_BIG_NUMBER;
      vU_b(5*k + 4) = update->gait[i*4 + j] * setup->f_max;
      k++;
    }
  }

  fpt mu = 1.f/setup->mu;
  Matrix<fpt,5,3> f_block;

  f_block <<  mu, 0,  1.f,
          -mu, 0,  1.f,
          0,  mu, 1.f,
          0, -mu, 1.f,
          0,   0, 1.f;

  for(s16 i = 0; i < setup->horizon*4; i++)
  {
    v_fmat.block(i*5,i*3,5,3) = f_block;
  }

  v_qH = 2*(vB_qp.transpose()*vS*vB_qp + update->alpha*v_eye_12h);
  v_qg = 2*vB_qp.transpose()*vS*(vA_qp*v_x_0 - vX_d);


  v_matrix_to_real(vH_qpoases,v_qH,setup->horizon*12, setup->horizon*12);
  v_matrix_to_real(vg_qpoases,v_qg,setup->horizon*12, 1);
  v_matrix_to_real(vA_qpoases,v_fmat,setup->horizon*20, setup->horizon*12);
  v_matrix_to_real(vub_qpoases,vU_b,setup->horizon*20, 1);

  for(s16 i = 0; i < 20*setup->horizon; i++)
    vlb_qpoases[i] = 0.0f;

  s16 num_constraints = 20*setup->horizon;
  s16 num_variables = 12*setup->horizon;


  qpOASES::int_t nWSR = 100;


  int new_vars = num_variables;
  int new_cons = num_constraints;

  for(int i =0; i < num_constraints; i++)
    v_con_elim[i] = 0;

  for(int i = 0; i < num_variables; i++)
    v_var_elim[i] = 0;


  for(int i = 0; i < num_constraints; i++)
  {
    if(! (v_near_zero(vlb_qpoases[i]) && v_near_zero(vub_qpoases[i]))) continue;
    double* c_row = &vA_qpoases[i*num_variables];
    for(int j = 0; j < num_variables; j++)
    {
      if(v_near_one(c_row[j]))
      {
        new_vars -= 3;
        new_cons -= 5;
        int cs = (j*5)/3 -3;
        v_var_elim[j-2] = 1;
        v_var_elim[j-1] = 1;
        v_var_elim[j  ] = 1;
        v_con_elim[cs] = 1;
        v_con_elim[cs+1] = 1;
        v_con_elim[cs+2] = 1;
        v_con_elim[cs+3] = 1;
        v_con_elim[cs+4] = 1;
      }
    }
  }
  //if(new_vars != num_variables)
  if(1==1)
  {
    int var_ind[new_vars];
    int v_con_ind[new_cons];
    int vc = 0;
    for(int i = 0; i < num_variables; i++)
    {
      if(!v_var_elim[i])
      {
        if(!(vc<new_vars))
        {
          printf("BAD ERROR 1\n");
        }
        var_ind[vc] = i;
        vc++;
      }
    }
    vc = 0;
    for(int i = 0; i < num_constraints; i++)
    {
      if(!v_con_elim[i])
      {
        if(!(vc<new_cons))
        {
          printf("BAD ERROR 1\n");
        }
        v_con_ind[vc] = i;
        vc++;
      }
    }
    for(int i = 0; i < new_vars; i++)
    {
      int olda = var_ind[i];
      vg_red[i] = vg_qpoases[olda];
      for(int j = 0; j < new_vars; j++)
      {
        int oldb = var_ind[j];
        vH_red[i*new_vars + j] = vH_qpoases[olda*num_variables + oldb];
      }
    }

    for (int con = 0; con < new_cons; con++)
    {
      for(int st = 0; st < new_vars; st++)
      {
        float cval = vA_qpoases[(num_variables*v_con_ind[con]) + var_ind[st] ];
        vA_red[con*new_vars + st] = cval;
      }
    }
    for(int i = 0; i < new_cons; i++)
    {
      int old = v_con_ind[i];
      vub_red[i] = vub_qpoases[old];
      vlb_red[i] = vlb_qpoases[old];
    }

    qpOASES::QProblem problem_red (new_vars, new_cons);
    qpOASES::Options op;
    op.setToMPC();
    op.printLevel = qpOASES::PL_NONE;
    problem_red.setOptions(op);
    //int_t nWSR = 50000;


    int rval = problem_red.init(vH_red, vg_red, vA_red, NULL, NULL, vlb_red, vub_red, nWSR);
    (void)rval;
    int rval2 = problem_red.getPrimalSolution(vq_red);
    if(rval2 != qpOASES::SUCCESSFUL_RETURN)
      printf("failed to solve!\n");

    // printf("solve time: %.3f ms, size %d, %d\n", solve_timer.getMs(), new_vars, new_cons);


    vc = 0;
    for(int i = 0; i < num_variables; i++)
    {
      if(v_var_elim[i])
      {
        vq_soln[i] = 0.0f;
      }
      else
      {
        vq_soln[i] = vq_red[vc];
        vc++;
      }
    }
  }
}
