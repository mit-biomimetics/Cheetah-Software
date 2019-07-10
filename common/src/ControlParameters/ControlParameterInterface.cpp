#include "ControlParameters/ControlParameterInterface.h"

std::string controlParameterRequestKindToString(
    ControlParameterRequestKind request) {
  switch (request) {
    case ControlParameterRequestKind::GET_USER_PARAM_BY_NAME:
      return "get user";
    case ControlParameterRequestKind::SET_USER_PARAM_BY_NAME:
      return "set user";
    case ControlParameterRequestKind::GET_ROBOT_PARAM_BY_NAME:
      return "get robot";
    case ControlParameterRequestKind::SET_ROBOT_PARAM_BY_NAME:
      return "set robot";
    default:
      return "unknown";
  }
}