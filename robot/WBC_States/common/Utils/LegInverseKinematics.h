#ifndef SPLINE_OPTIMIZER
#define SPLINE_OPTIMIZER

#include <nlopt.hpp>
#include <math.h>
#include "Utilities/BSplineBasic.h"
#include "Utilities/Utilities_print.h"

// Initial: (pos, vel, acc) constrained
// Final: (pos, vel) constrained, acceleration is free
//
template<typename T, int DIM, int NUM_MIDDLE>
class SplineOptimizer{
    public:
        SplineOptimizer():_num_check_pt(4){
            //_mid_pt = new T*[NUM_MIDDLE];
            for(size_t i(0); i<NUM_MIDDLE; ++i){
                _mid_pt[i] = new T[DIM];
            }
        }
        ~SplineOptimizer(){
            for(size_t i(0); i<NUM_MIDDLE; ++i){
                delete [] _mid_pt[i];
            }
        }
        static double MinAcc(
                const std::vector<double> &x, 
                std::vector<double> & grad,
                void *d_) {

            (void)grad;
            (void)d_;

            T cost_sum(0.);
            for(size_t i(0); i<DIM; ++i){
                //grad[i] = 2.*x[i]; 
                //cost_sum += (x[i]*x[i]); 
                cost_sum += (x[i]); 
            }
            return cost_sum;
        } 

        static void BoundAcceleration(
                unsigned m, double * result, unsigned n, const double *x, 
                double * grad, void*data){

            (void)m;
            (void)n;
            (void)grad;
            //for(size_t i(0); i< m; ++i) grad[i] = 0.;

            SplineOptimizer* tester = (SplineOptimizer*)data;
            T delta_t = tester->_fin_time/((T)tester->_num_check_pt);

            for (unsigned int k(0); k<DIM; ++k) tester->_fin_opt[DIM*2 + k] = x[k + DIM];

            tester->_spline.SetParam(tester->_ini, tester->_fin_opt, 
                    tester->_mid_pt, tester->_fin_time);

            //printf("ini pos: %f, %f\n", tester->_ini[0], tester->_ini[1]);
            //printf("ini vel: %f, %f\n", tester->_ini[2], tester->_ini[3]);
            //printf("ini acc: %f, %f\n", tester->_ini[4], tester->_ini[5]);

            //printf("fin pos: %f, %f\n", tester->_fin_opt[0], tester->_fin_opt[1]);
            //printf("fin vel: %f, %f\n", tester->_fin_opt[2], tester->_fin_opt[3]);
            //printf("fin acc: %f, %f\n", tester->_fin_opt[4], tester->_fin_opt[5]);

            //printf("mid pt: %f, %f\n", tester->_mid_pt[0][0], tester->_mid_pt[0][1]);


            T acc[DIM]; 
            for(size_t i(0); i< DIM; ++i) acc[i] = 0.;

            for(unsigned int i(0); i<tester->_num_check_pt; ++i){
                tester->_spline.getCurveDerPoint(delta_t * (T)(i+1), 2, acc);
                //printf("%f sec acc: %f, %f\n", delta_t*(i+1), acc[0], acc[1]);
                for(int j(0); j<DIM; ++j){
                    result[2*DIM*i + 2*j] = acc[j] - x[j];
                    result[2*DIM*i + 2*j+1] = -acc[j] - x[j];
                }
            }
            //pretty_print(result, "min_max_acc", m);
            //pretty_print(x, "opt var", n);
        }

        void setParam(T* ini, T* fin, T ** middle_pt, T fin_time){
            for(size_t i(0); i<DIM*3; ++i) _fin_opt[i] = fin[i];
            for(size_t i(0); i<DIM*3; ++i) _ini[i]= ini[i];
            for(size_t i(0); i<NUM_MIDDLE; ++i) 
                for(size_t j(0); j<DIM; ++j)
                    _mid_pt[i][j]= middle_pt[i][j];

            _fin_time = fin_time;

            // (bound X DIM), (Final Acceleration X DIM)
            int num_opt_var = 2*DIM;
            nlopt::opt* test = new nlopt::opt(nlopt::LN_COBYLA, num_opt_var);
            //nlopt::opt* test = new nlopt::opt(nlopt::LD_MMA, num_opt_var);
            //nlopt::opt* test = new nlopt::opt(nlopt::LD_LBFGS, num_opt_var);
            //nlopt::opt* test = new nlopt::opt(nlopt::GD_STOGO, num_opt_var);
            //nlopt::opt* test = new nlopt::opt(nlopt::LD_SLSQP, num_opt_var);
            //nlopt::opt* test = new nlopt::opt(nlopt::LD_CCSAQ, num_opt_var);
            std::vector<double> x(num_opt_var);

            std::vector<double> lb(num_opt_var);
            std::vector<double> ub(num_opt_var);
            for(int i(0);i<DIM; ++i){
                lb[i] = 0.;
                ub[i] = HUGE_VAL;
            }
            for(int i(0);i<DIM; ++i){
                lb[i + DIM] = -HUGE_VAL;
                ub[i + DIM] = HUGE_VAL;
            }
            test->set_lower_bounds(lb);
            test->set_upper_bounds(ub);
            
            int num_acc_check_pt = 2 * (DIM*_num_check_pt);
            std::vector<double> tol_ieq(num_acc_check_pt);
            for(int i(0); i<num_acc_check_pt; ++i){ tol_ieq[i] = 1e-4; }
            test->add_inequality_mconstraint(SplineOptimizer::BoundAcceleration, this, tol_ieq);

            // Objective Function
            test->set_xtol_rel(1e-10);
            test->set_ftol_rel(1e-10);
            test->set_maxeval(500); 

            test->set_min_objective(SplineOptimizer::MinAcc, this);
            T opt_f;
            nlopt::result opt_result = test->optimize(x, opt_f);

            //printf("result, cost: %d, %f\n", opt_result, opt_f);
            //pretty_print(x, "opt result");
            (void)opt_result;

            //std::cout<<opt_result<<std::endl;

            for(size_t i(0); i<DIM; ++i) _fin_opt[i + DIM*2] = x[i+ DIM];
            _spline.SetParam(ini, _fin_opt, middle_pt, fin_time);
        }
        void getPos(T time, T* pos){
            _spline.getCurvePoint(time, pos);
        }
        void getVel(T time, T* vel){
            _spline.getCurveDerPoint(time, 1, vel);
        }
        void getAcc(T time, T* acc){
            _spline.getCurveDerPoint(time, 2, acc);
        }



        size_t _num_check_pt;
        T _fin_time;
        T _fin_opt[DIM*3];
        T _ini[DIM*3];
        //T** _mid_pt;
        T* _mid_pt[NUM_MIDDLE];
        BS_Basic<T, DIM, 3, NUM_MIDDLE, 2, 2> _spline;
};
#endif
