#ifndef FILTERS_
#define FILTERS_

template <typename T>
class filter {
 public:
  filter(void) {}
  virtual ~filter(void) {}
  virtual void input(T input_value) = 0;
  virtual T output(void) = 0;
  virtual void clear(void) = 0;
};

template <typename T>
class butterworth_filter : public filter<T> {
 public:
  butterworth_filter(int num_sample, T dt, T cutoff_frequency);
  virtual ~butterworth_filter(void);
  virtual void input(T input_value);
  virtual T output(void);
  virtual void clear(void);

 private:
  T *mpBuffer;
  int mCurIdx;
  int mNumSample;
  T mDt;
  T mCutoffFreq;
  T mValue;
};

template <typename T>
class digital_lp_filter : public filter<T> {
 public:
  digital_lp_filter(T w_c, T t_s);
  virtual ~digital_lp_filter(void);
  virtual void input(T input_value);
  virtual T output(void);
  virtual void clear(void);

 private:
  T Lpf_in_prev[2];
  T Lpf_out_prev[2];
  T Lpf_in1, Lpf_in2, Lpf_in3, Lpf_out1, Lpf_out2;
  T lpf_out;
};

template <typename T>
class moving_average_filter : public filter<T> {
 public:
  moving_average_filter(int num_data);
  virtual ~moving_average_filter();
  virtual void input(T input_value);
  virtual T output(void);
  virtual void clear(void);

 private:
  T *buffer_;
  int num_data_;
  int idx_;
  T sum_;
};

template <typename T>
class deriv_lp_filter : public filter<T> {
 public:
  deriv_lp_filter(T w_c, T t_s);
  virtual ~deriv_lp_filter(void);
  virtual void input(T input_value);
  virtual T output(void);
  virtual void clear(void);

 private:
  T Lpf_in_prev[2];
  T Lpf_out_prev[2];
  T Lpf_in1, Lpf_in2, Lpf_in3, Lpf_out1, Lpf_out2;
  T lpf_out;
};

template <typename T>
class ff01_filter : public filter<T> {
 public:
  ff01_filter(float t_s, float w_c);
  virtual ~ff01_filter(void);
  virtual void input(T input_value);
  virtual T output(void);
  virtual void clear(void);

 private:
  T Lpf_in_prev[2];
  T Lpf_out_prev[2];
  T Lpf_in1, Lpf_in2, Lpf_in3, Lpf_out1, Lpf_out2;
  T lpf_out;
};

template <typename T>
class ff02_filter : public filter<T> {
 public:
  ff02_filter(float t_s, float w_c);
  virtual ~ff02_filter(void);
  virtual void input(T input_value);
  virtual T output(void);
  virtual void clear(void);

 private:
  T Lpf_in_prev[2];
  T Lpf_out_prev[2];
  T Lpf_in1, Lpf_in2, Lpf_in3, Lpf_out1, Lpf_out2;
  T lpf_out;
};

template <typename T>
class AverageFilter : public filter<T> {
 public:
  AverageFilter(T dt, T t_const, T limit);
  virtual ~AverageFilter();
  virtual void input(T input_value);
  virtual T output(void);
  virtual void clear(void);

 private:
  T est_value_;
  T dt_;
  T t_const_;
  T limit_;
};
#endif
