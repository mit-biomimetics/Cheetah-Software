/*! @file ControlParameterInterface.cpp
 *  @brief Types to allow remote access of control parameters, for use with
 * LCM/Shared memory
 *
 * There are response/request messages.  The robot receives request messages and
 * responds with response messages The request messages either request setting a
 * parameter or getting a parameter
 */

#include "ControlParameters/ControlParameterInterface.h"

/*!
 * Convert control parameter request kind into a string
 */
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