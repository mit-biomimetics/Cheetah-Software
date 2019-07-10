#ifndef PATH_H
#define PATH_H

#include <Utilities/Utilities_print.h>
#include <cppTypes.h>

template <typename T>
class Path {
 public:
  Path() {}
  virtual ~Path() {}

  virtual void getDesLoc(const T& t, const Vec3<T>& curr_loc,
                         Vec3<T>& des_loc) = 0;
};

template <typename T>
class LinPath : public Path<T> {
 public:
  LinPath()
      : _b_first_visit(true), _amp(0.6), _freq(0.06), _prepare_time(2.0) {}
  virtual ~LinPath() {}

  virtual void getDesLoc(const T& t, const Vec3<T>& curr_loc,
                         Vec3<T>& des_loc) {
    if (_b_first_visit) {
      _ini_time = t;
      _b_first_visit = false;
    }
    T path_time = t - _ini_time;
    // printf("time: %f\n", path_time);
    if (path_time < _prepare_time) {
      _target_loc = curr_loc;
      _ini_loc = curr_loc;
    } else {
      _target_loc[0] = _ini_loc[0] + _amp * sin(2. * M_PI * _freq *
                                                (path_time - _prepare_time));
    }
    des_loc = _target_loc;
  }

 protected:
  T _ini_time;
  bool _b_first_visit;
  T _amp, _freq, _prepare_time;
  Vec3<T> _target_loc, _ini_loc;
};

template <typename T>
class CircularPath : public Path<T> {
 public:
  CircularPath()
      : _loop_time(0.),
        _prepare_time(1.0),
        _b_first_visit(true),
        _radius_x(0.9),
        _radius_y(0.25),
        _freq(0.08) {}
  virtual ~CircularPath() {}

  virtual void getDesLoc(const T& t, const Vec3<T>& curr_loc,
                         Vec3<T>& des_loc) {
    if (_b_first_visit) {
      _ini_time = t;
      _b_first_visit = false;
    }
    T path_time = t - _ini_time;

    if (path_time < _prepare_time) {
      _target_loc = curr_loc;
      _ini_loc = curr_loc;
    } else {
      _target_loc[0] =
          _ini_loc[0] +
          _radius_x * sin(2. * M_PI * _freq * (path_time - _prepare_time));
      _target_loc[1] =
          _ini_loc[1] +
          _radius_y *
              (cos(2. * M_PI * _freq * (path_time - _prepare_time)) - 1.);

      _target_loc[2] = -2. * M_PI * (path_time - _prepare_time) * _freq;
      // TEST
      //_target_loc[2] = 0.;
    }
    des_loc = _target_loc;
  }

 protected:
  T _ini_time, _loop_time, _prepare_time;
  bool _b_first_visit;
  T _radius_x, _radius_y, _freq;
  Vec3<T> _target_loc, _ini_loc;
};

#endif
