#ifndef PROJECT_UTILITIES_H
#define PROJECT_UTILITIES_H

#include <algorithm>
#include <map>
#include <random>
#include <unordered_map>
#include <vector>
#include "cppTypes.h"

/*!
 * Are two floating point values almost equal?
 */
template <typename T>
bool fpEqual(T a, T b, T tol) {
  return std::abs(a - b) <= tol;
}

/*!
 * Are two std::vectors equal?
 */
template <typename T>
bool vectorEqual(const std::vector<T>& a, const std::vector<T>& b) {
  if (a.size() != b.size()) return false;
  for (size_t i = 0; i < a.size(); i++) {
    if (a[i] != b[i]) return false;
  }
  return true;
}

/*!
 * Coerce in to be between min and max
 */
template <typename T>
T coerce(T in, T min, T max) {
  if (in < min) {
    in = min;
  }
  if (in > max) {
    in = max;
  }
  return in;
}

template <typename T>
T deadband(T x, T range) {
  if (x < range && x > -range) x = T(0);
  return x;
}

template <typename T>
void eigenDeadband(Eigen::MatrixBase<T>& v, typename T::Scalar band) {
  for (size_t i = 0; i < T::RowsAtCompileTime; i++) {
    for (size_t j = 0; j < T::ColsAtCompileTime; j++) {
      v(i, j) = deadband(v(i, j), band);
    }
  }
}

/*!
 * Get the sign of a number
 * 1 for positive, 0 for 0, -1 for negative...
 */
template <typename T>
int sgn(T val) {
  return (T(0) < val) - (val < T(0));
}

/*!
 * Fill an eigen type with random numbers from a random generator and uniform
 * real distribution.
 * TODO: is there a way to make this work nicely with normal distributions too?
 */
template <typename T>
void fillEigenWithRandom(
    Eigen::MatrixBase<T>& v, std::mt19937& gen,
    std::uniform_real_distribution<typename T::Scalar>& dist) {
  for (size_t i = 0; i < T::RowsAtCompileTime; i++) {
    for (size_t j = 0; j < T::ColsAtCompileTime; j++) {
      v(i, j) = dist(gen);
    }
  }
}

/*!
 * Generate a random number following normal distribution
 */
template <typename T>
T generator_gaussian_noise(T mean, T var) {
  static bool hasSpare = false;
  static T rand1, rand2;

  if (hasSpare) {
    hasSpare = false;
    return mean + sqrt(var * rand1) * sin(rand2);
  }
  hasSpare = true;

  rand1 = rand() / ((T)RAND_MAX);
  if (rand1 < 1e-100) rand1 = 1e-100;
  rand1 = -2 * log(rand1);
  rand2 = rand() / ((T)RAND_MAX) * M_PI * 2.;

  // printf("rand: %f, %f\n", rand1, rand2);
  return mean + sqrt(var * rand1) * cos(rand2);
}

/*!
 * Does the unordered map contain the given element?
 */
template <typename T1, typename T2>
bool uMapContains(const std::unordered_map<T1, T2>& set, T1 key) {
  return set.find(key) != set.end();
}

/*!
 * Does the unordered map contain the given element?
 */
template <typename T1, typename T2>
bool mapContains(const std::map<T1, T2>& set, T1 key) {
  return set.find(key) != set.end();
}

/*!
 * Convert a floating point number to a string.  Is preferable over
 * std::to_string because this uses scientific notation and won't truncate
 * small/large numbers.
 */
template <typename T>
std::string numberToString(T number) {
  static_assert(std::is_floating_point<T>::value,
                "numberToString must use a floating point type!");
  char buffer[100];
  sprintf(buffer, "%g", number);
  return std::string(buffer);
}

/*!
 * map value x in (inputMin, inputMax) to (outputMin, outputMax)
 */
template <typename T>
T mapToRange(T x, T inputMin, T inputMax, T outputMin, T outputMax) {
  return outputMin +
         (x - inputMin) * (outputMax - outputMin) / (inputMax - inputMin);
}

template <typename T>
std::string eigenToString(Eigen::MatrixBase<T>& value) {
  std::stringstream ss;
  ss << value;
  return ss.str();
}

static inline std::string boolToString(bool b) {
  return std::string(b ? "true" : "false");
}

void writeStringToFile(const std::string& fileName,
                       const std::string& fileData);
std::string getCurrentTimeAndDate();
std::string getConfigDirectoryPath();

/*!
 * Get the rotation matrix coincide with euler angle
 * Intrisic ZYX rotation
 */
template <typename T>
void EulerZYX_2_SO3(const Vec3<T>& euler_zyx, Mat3<T>& SO3) {
  Mat3<T> Mat3_Z, Mat3_Y, Mat3_X;
  Mat3_Z << cos(euler_zyx[0]), -sin(euler_zyx[0]), 0, sin(euler_zyx[0]),
      cos(euler_zyx[0]), 0, 0, 0, 1;
  Mat3_Y << cos(euler_zyx[1]), 0, sin(euler_zyx[1]), 0, 1, 0,
      -sin(euler_zyx[1]), 0, cos(euler_zyx[1]);
  Mat3_X << 1, 0, 0, 0, cos(euler_zyx[2]), -sin(euler_zyx[2]), 0,
      sin(euler_zyx[2]), cos(euler_zyx[2]);

  SO3 = Mat3_Z * Mat3_Y * Mat3_X;
}

// Smooth Changing
template <typename T>
T smooth_change(T ini, T end, T moving_duration, T curr_time) {
  if (curr_time > moving_duration) {
    return end;
  }
  return (ini +
          (end - ini) * 0.5 * (1 - cos(curr_time / moving_duration * M_PI)));
}

template <typename T>
T smooth_change_vel(T ini, T end, T moving_duration, T curr_time) {
  if (curr_time > moving_duration) {
    return 0.0;
  }
  return ((end - ini) * 0.5 * (M_PI / moving_duration) *
          sin(curr_time / moving_duration * M_PI));
}

template <typename T>
T smooth_change_acc(T ini, T end, T moving_duration, T curr_time) {
  if (curr_time > moving_duration) {
    return 0.0;
  }
  return ((end - ini) * 0.5 * (M_PI / moving_duration) *
          (M_PI / moving_duration) * cos(curr_time / moving_duration * M_PI));
}

template <typename T>
T stringToNumber(const std::string& str) {
  static_assert(std::is_same<T, double>::value || std::is_same<T, float>::value,
                "stringToNumber only works for double/float");

  if (std::is_same<T, double>::value) {
    return std::stod(str);
  } else if (std::is_same<T, float>::value) {
    return std::stof(str);
  }
}

template <typename T>
T stringToNumber(const char* str) {
  return stringToNumber<T>(std::string(str));
}

template <typename T>
Vec3<T> stringToVec3(const std::string& str) {
  Vec3<T> v;
  size_t i = 0;

  // seek past whitespace
  while (str.at(i) == ' ') i++;

  if (str.at(i) == '[') {
    i++;
  } else {
    throw std::runtime_error("stringToVec3 didn't find open bracket");
  }

  // seek past whitespace
  while (str.at(i) == ' ') i++;
  size_t start = i;

  // seek to end of first number
  while (str.at(i) != ',') i++;
  v[0] = stringToNumber<T>(str.substr(start, i - start));
  i++;

  while (str.at(i) == ' ') i++;
  start = i;
  while (str.at(i) != ',') i++;
  v[1] = stringToNumber<T>(str.substr(start, i - start));
  i++;

  while (str.at(i) == ' ') i++;
  start = i;
  while (str.at(i) != ']') i++;
  v[2] = stringToNumber<T>(str.substr(start, i - start));
  return v;
}

std::string getLcmUrl(s64 ttl);

#endif  // PROJECT_UTILITIES_H
