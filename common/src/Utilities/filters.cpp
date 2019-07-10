#include "include/Utilities/filters.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define BUDDA_Q_SCALE 6.f

template <typename T>
moving_average_filter<T>::moving_average_filter(int num_data)
    : num_data_(num_data), idx_(0), sum_(0.0) {
  buffer_ = new T[num_data_];
  memset((void *)buffer_, 0.0, sizeof(T) * num_data_);
}

template <typename T>
void moving_average_filter<T>::input(T input_value) {
  sum_ -= buffer_[idx_];
  sum_ += input_value;
  buffer_[idx_] = input_value;
  ++idx_;
  idx_ %= num_data_;
}

template <typename T>
T moving_average_filter<T>::output() {
  return sum_ / num_data_;
}

template <typename T>
void moving_average_filter<T>::clear(void) {
  sum_ = 0.0;
  memset((void *)buffer_, 0.0, sizeof(T) * num_data_);
}

template <typename T>
moving_average_filter<T>::~moving_average_filter() {
  delete[] buffer_;
}

template class moving_average_filter<double>;
template class moving_average_filter<float>;

/*============================================================================*/

template <typename T>
butterworth_filter<T>::butterworth_filter(int num_sample, T dt,
                                          T cutoff_frequency) {
  mNumSample = num_sample;
  mDt = dt;
  mCutoffFreq = cutoff_frequency;

  mpBuffer = new T[num_sample];
  memset((void *)mpBuffer, 0, sizeof(T) * num_sample);

  mCurIdx = 0;
}

template <typename T>
butterworth_filter<T>::~butterworth_filter(void) {
  delete[] mpBuffer;
}

template <typename T>
void butterworth_filter<T>::input(T input_value) {
  int j;
  T sqrt_2 = sqrt(2);
  T value = 0;
  for (j = mNumSample - 2; j >= 0; j--) {
    mpBuffer[j + 1] = mpBuffer[j];
  }

  mpBuffer[0] = input_value;
  for (j = 0; j < mNumSample; j++) {
    T t = (T)j * mDt;
    value += sqrt_2 / mCutoffFreq * mpBuffer[j] * exp(-1. / sqrt_2 * t) *
             sin(mCutoffFreq / sqrt_2 * t) * mDt;
    //		value += sqrt_2 * exp(-1./sqrt_2*t) * sin(1./sqrt_2*t ) * mDt;
  }
  mValue = value;
}

template <typename T>
T butterworth_filter<T>::output(void) {
  return mValue;
}

template <typename T>
void butterworth_filter<T>::clear(void) {
  for (int i(0); i < mNumSample; ++i) {
    mpBuffer[i] = 0.0;
  }
}

template class butterworth_filter<double>;
template class butterworth_filter<float>;

/*============================================================================*/

template <typename T>
digital_lp_filter<T>::digital_lp_filter(T w_c, T t_s) {
  Lpf_in_prev[0] = Lpf_in_prev[1] = 0;
  Lpf_out_prev[0] = Lpf_out_prev[1] = 0;
  Lpf_in1 = 0, Lpf_in2 = 0, Lpf_in3 = 0, Lpf_out1 = 0, Lpf_out2 = 0;
  float den = 2500 * t_s * t_s * w_c * w_c + 7071 * t_s * w_c + 10000;

  Lpf_in1 = 2500 * t_s * t_s * w_c * w_c / den;
  Lpf_in2 = 5000 * t_s * t_s * w_c * w_c / den;
  Lpf_in3 = 2500 * t_s * t_s * w_c * w_c / den;
  Lpf_out1 = -(5000 * t_s * t_s * w_c * w_c - 20000) / den;
  Lpf_out2 = -(2500 * t_s * t_s * w_c * w_c - 7071 * t_s * w_c + 10000) / den;
}

template <typename T>
digital_lp_filter<T>::~digital_lp_filter(void) {}

template <typename T>
void digital_lp_filter<T>::input(T lpf_in) {
  lpf_out = Lpf_in1 * lpf_in + Lpf_in2 * Lpf_in_prev[0] +
            Lpf_in3 * Lpf_in_prev[1] +  // input component
            Lpf_out1 * Lpf_out_prev[0] +
            Lpf_out2 * Lpf_out_prev[1];  // output component
  Lpf_in_prev[1] = Lpf_in_prev[0];
  Lpf_in_prev[0] = lpf_in;
  Lpf_out_prev[1] = Lpf_out_prev[0];
  Lpf_out_prev[0] = lpf_out;
}

template <typename T>
T digital_lp_filter<T>::output(void) {
  return lpf_out;
}

template <typename T>
void digital_lp_filter<T>::clear(void) {
  Lpf_in_prev[1] = 0;
  Lpf_in_prev[0] = 0;
  Lpf_out_prev[1] = 0;
  Lpf_out_prev[0] = 0;
}

template class digital_lp_filter<double>;
template class digital_lp_filter<float>;

/*============================================================================*/

template <typename T>
deriv_lp_filter<T>::deriv_lp_filter(T w_c, T t_s) {
  Lpf_in_prev[0] = 0;
  Lpf_in_prev[1] = 0;
  Lpf_out_prev[0] = 0;
  Lpf_out_prev[1] = 0;
  Lpf_in1 = 0;
  Lpf_in2 = 0;
  Lpf_in3 = 0;
  Lpf_out1 = 0;
  Lpf_out2 = 0;
  T a = 1.4142;
  T den = 4 + 2 * a * w_c * t_s + t_s * t_s * w_c * w_c;

  Lpf_in1 = 2 * t_s * w_c * w_c / den;
  Lpf_in2 = 0;
  Lpf_in3 = -2. * t_s * w_c * w_c / den;
  Lpf_out1 = -1. * (-8 + t_s * t_s * w_c * w_c * 2) / den;
  Lpf_out2 = -1. * (4 - 2 * a * w_c * t_s + t_s * t_s * w_c * w_c) / den;
  lpf_out = 0.0;
  clear();
}

template <typename T>
deriv_lp_filter<T>::~deriv_lp_filter(void) {}

template <typename T>
void deriv_lp_filter<T>::input(T lpf_in) {
  // static int i(0);
  lpf_out = Lpf_in1 * lpf_in + Lpf_in2 * Lpf_in_prev[0] +
            Lpf_in3 * Lpf_in_prev[1] +  // input component
            Lpf_out1 * Lpf_out_prev[0] +
            Lpf_out2 * Lpf_out_prev[1];  // output component

  // printf("%i th filter (%f): %f, %f, %f, %f, %f, %f, %f, %f, %f, %f\n",i,
  //        lpf_out,
  //        Lpf_in1, lpf_in, Lpf_in2,
  //        Lpf_in_prev[0], Lpf_in3, Lpf_in_prev[1],
  //        Lpf_out1, Lpf_out_prev[0], Lpf_out2, Lpf_out_prev[1]);

  // if(lpf_out>100){
  //     exit(0);
  // }

  Lpf_in_prev[1] = Lpf_in_prev[0];
  Lpf_in_prev[0] = lpf_in;
  Lpf_out_prev[1] = Lpf_out_prev[0];
  Lpf_out_prev[0] = lpf_out;
  // ++i;
}

template <typename T>
T deriv_lp_filter<T>::output(void) {
  return lpf_out;
}

template <typename T>
void deriv_lp_filter<T>::clear(void) {
  Lpf_in_prev[1] = 0;
  Lpf_in_prev[0] = 0;
  Lpf_out_prev[1] = 0;
  Lpf_out_prev[0] = 0;
}

template class deriv_lp_filter<double>;
template class deriv_lp_filter<float>;

/*============================================================================*/

template <typename T>
ff01_filter<T>::ff01_filter(float t_s, float w_c) {
  Lpf_in_prev[0] = Lpf_in_prev[1] = 0;
  Lpf_out_prev[0] = Lpf_out_prev[1] = 0;
  Lpf_in1 = 0, Lpf_in2 = 0, Lpf_in3 = 0, Lpf_out1 = 0, Lpf_out2 = 0;
  T a = 1.4142;
  T den = 4 + 2 * a * w_c * t_s + t_s * t_s * w_c * w_c;
  T J = 0.00008;
  T B = 0.0002;

  Lpf_in1 = B * t_s * t_s * w_c * w_c + 2 * J * t_s * w_c * w_c;
  Lpf_in2 = 2 * B * t_s * t_s * w_c * w_c;
  Lpf_in3 = B * t_s * t_s * w_c * w_c - 2 * J * t_s * w_c * w_c;
  Lpf_out1 = -1. * (-8 + t_s * t_s * w_c * w_c * 2) / den;
  Lpf_out2 = -1. * (4 - 2 * a * w_c * t_s + t_s * t_s * w_c * w_c) / den;
}

template <typename T>
ff01_filter<T>::~ff01_filter(void) {}

template <typename T>
void ff01_filter<T>::input(T lpf_in) {
  lpf_out = Lpf_in1 * lpf_in + Lpf_in2 * Lpf_in_prev[0] +
            Lpf_in3 * Lpf_in_prev[1] +  // input component
            Lpf_out1 * Lpf_out_prev[0] +
            Lpf_out2 * Lpf_out_prev[1];  // output component
  Lpf_in_prev[1] = Lpf_in_prev[0];
  Lpf_in_prev[0] = lpf_in;
  Lpf_out_prev[1] = Lpf_out_prev[0];
  Lpf_out_prev[0] = lpf_out;
}

template <typename T>
T ff01_filter<T>::output(void) {
  return lpf_out;
}

template <typename T>
void ff01_filter<T>::clear(void) {
  Lpf_in_prev[1] = 0;
  Lpf_in_prev[0] = 0;
  Lpf_out_prev[1] = 0;
  Lpf_out_prev[0] = 0;
}

template class ff01_filter<float>;
template class ff01_filter<double>;

/*============================================================================*/

template <typename T>
ff02_filter<T>::ff02_filter(float t_s, float w_c) {
  T J = 0.003216;

  Lpf_in_prev[0] = Lpf_in_prev[1] = 0;
  Lpf_out_prev[0] = Lpf_out_prev[1] = 0;
  Lpf_in1 = 0, Lpf_in2 = 0, Lpf_in3 = 0, Lpf_out1 = 0, Lpf_out2 = 0;

  T a = 1.4142;
  T den = 4 + 2 * a * w_c * t_s + t_s * t_s * w_c * w_c;

  Lpf_in1 = J * 2 * t_s * w_c * w_c / den;
  Lpf_in2 = 0;
  Lpf_in3 = -2. * J * t_s * w_c * w_c / den;
  Lpf_out1 = -1. * (-8 + t_s * t_s * w_c * w_c * 2) / den;
  Lpf_out2 = -1. * (4 - 2 * a * w_c * t_s + t_s * t_s * w_c * w_c) / den;

  clear();
}

template <typename T>
ff02_filter<T>::~ff02_filter(void) {}

template <typename T>
void ff02_filter<T>::input(T lpf_in) {
  lpf_out = Lpf_in1 * lpf_in + Lpf_in2 * Lpf_in_prev[0] +
            Lpf_in3 * Lpf_in_prev[1] +  // input component
            Lpf_out1 * Lpf_out_prev[0] +
            Lpf_out2 * Lpf_out_prev[1];  // output component
  Lpf_in_prev[0] = lpf_in;
  Lpf_in_prev[1] = Lpf_in_prev[0];
  Lpf_out_prev[0] = lpf_out;
  Lpf_out_prev[1] = Lpf_out_prev[0];
}

template <typename T>
T ff02_filter<T>::output(void) {
  return lpf_out;
}

template <typename T>
void ff02_filter<T>::clear(void) {
  Lpf_in_prev[1] = 0;
  Lpf_in_prev[0] = 0;
  Lpf_out_prev[1] = 0;
  Lpf_out_prev[0] = 0;
}

template class ff02_filter<float>;
template class ff02_filter<double>;

/*============================================================================*/

template <typename T>
AverageFilter<T>::AverageFilter(T dt, T t_const, T limit)
    : dt_(dt), t_const_(t_const), limit_(limit) {
  est_value_ = 0.;
}

template <typename T>
AverageFilter<T>::~AverageFilter() {
  est_value_ = 0;
}

template <typename T>
void AverageFilter<T>::clear() {
  est_value_ = 0.;
}

template <typename T>
void AverageFilter<T>::input(T input) {
  T update_value = input - est_value_;
  if (fabs(update_value) > limit_) {
    update_value = 0.;
  }
  est_value_ += (dt_ / (dt_ + t_const_)) * update_value;
}

template <typename T>
T AverageFilter<T>::output() {
  return est_value_;
}

template class AverageFilter<float>;
template class AverageFilter<double>;
