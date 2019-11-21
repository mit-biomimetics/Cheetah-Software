#ifndef QPSOLVER_PROBLEMGENERATOR_H
#define QPSOLVER_PROBLEMGENERATOR_H

#include "QpProblem.h"

template<typename T>
class ProblemGenerator
{
public:
  QpProblem<T> generateSparseMPC(s64 n_states, s64 n_control, s64 horizon, s64 n_constraint);
};


#endif //QPSOLVER_PROBLEMGENERATOR_H
