/*! @file FirstOrderIIRFilter.h
 *  @brief A simple first order filter
 *
 *
 */

#ifndef PROJECT_FIRSTORDERIIRFILTER_H
#define PROJECT_FIRSTORDERIIRFILTER_H

#include <cmath>

/*!
 * First Order Filter
 * @tparam T : type of the data to be filtered
 * @tparam T2 : floating point type for the cutoff/sample frequencies, gain
 */
template <typename T, typename T2>
class FirstOrderIIRFilter {
 public:
  /*!
   * Create a new first order filter
   * @param cutoffFrequency : cutoff frequency of filter
   * @param sampleFrequency : sample frequency (rate at which update is called)
   * @param initialialValue : initial value stored in the filter
   */
  FirstOrderIIRFilter(T2 cutoffFrequency, T2 sampleFrequency, T& initialValue) {
    _alpha = 1 - std::exp(-2 * M_PI * cutoffFrequency / sampleFrequency);
    _state = initialValue;
  }

  /*!
   * Create a new first order filter
   * @param alpha : filter parameter
   * @param initialValue : initial value
   */
  FirstOrderIIRFilter(T2 alpha, T& initialValue)
      : _state(initialValue), _alpha(alpha) {}

  /*!
   * Update the filter with a new sample
   * @param x : the new sample
   * @return the new state of the filter
   */
  T update(T& x) {
    _state = _alpha * x + (T2(1) - _alpha) * _state;
    return _state;
  }

  /*!
   * Get the value of the filter, without updating
   */
  T get() { return _state; }

  void reset() { _state *= T2(0); }

 private:
  T _state;
  T2 _alpha;
};

#endif  // PROJECT_FIRSTORDERIIRFILTER_H
