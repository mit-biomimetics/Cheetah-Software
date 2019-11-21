#ifndef PROJECT_EDGETRIGGER_H
#define PROJECT_EDGETRIGGER_H

template<typename T>
class EdgeTrigger {
public:
  EdgeTrigger(T initial_state) : _state(initial_state) { }

  bool trigger(T& x) {
    if(_state == x) {
      return false;
    } else {
      _state = x;
      return true;
    }
  }

private:
  T _state;
};

#endif //PROJECT_EDGETRIGGER_H
