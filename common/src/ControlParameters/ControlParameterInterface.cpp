#include "ControlParameters/ControlParameterInterface.h"

std::string controlParameterRequestKindToString(ControlParameterRequestKind request) {
  switch(request) {
    case ControlParameterRequestKind::GET_PARAM_BY_NAME:
      return "get";
    case ControlParameterRequestKind::SET_PARAM_BY_NAME:
      return "set";
    default:
      return "unknown";
  }
}