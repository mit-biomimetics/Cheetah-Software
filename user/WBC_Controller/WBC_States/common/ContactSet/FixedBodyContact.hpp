#ifndef Cheetah_FIXED_BODY_CONTACT
#define Cheetah_FIXED_BODY_CONTACT

#include <WBC/ContactSpec.hpp>

template <typename T>
class FixedBodyContact : public ContactSpec<T> {
 public:
  FixedBodyContact();
  virtual ~FixedBodyContact();

 protected:
  virtual bool _UpdateJc();
  virtual bool _UpdateJcDotQdot();
  virtual bool _UpdateUf();
  virtual bool _UpdateInequalityVector();
};

#endif
