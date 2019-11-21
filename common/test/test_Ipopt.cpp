#include "gmock/gmock.h"
#include "gtest/gtest.h"

#ifdef IPOPT_OPTION 

#include <IpTNLP.hpp>
using namespace Ipopt;
class JumpNLP: public TNLP{
};

TEST(Ipopt, JumpNLP){
  printf("Jump NLP start\n");

}

#endif
