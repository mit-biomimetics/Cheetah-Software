/*! @file ControlParameters.h
 *  @brief Interface to set gains/control parameters for simulator and robot
 *  These are designed to be updated infrequently.  For high frequency data,
 * consider using Driver Inputs or adding to the Robot Debug Data instead.
 *
 * ControlParameter: a single value, either a double, float, or s64.  Each
 * control parameter must have a unique name Control parameters know their type
 * as well as if they've been initialized or not ControlParameters must be
 * initialized, either from reading from a file, reading from LCM, or some other
 * way (TODO consider supporting Vec3's?) All control parameters must go in a
 * ControlParameters (for now, just robot settings and sim settings)
 *
 * TODO - what should happen when this fails?
 *
 * See test_ControlParameters for an example of how this works
 */

#ifndef PROJECT_CONTROLPAREMETERS_H
#define PROJECT_CONTROLPAREMETERS_H

#include <map>
#include <mutex>
#include <string>
#include "Utilities/utilities.h"
#include "cTypes.h"

#define CONTROL_PARAMETER_MAXIMUM_NAME_LENGTH 64

#define INIT_PARAMETER(name) param_##name(#name, name, collection)
#define DECLARE_PARAMETER(type, name) \
  type name;                          \
  ControlParameter param_##name;

enum class ControlParameterValueKind : u64 {
  FLOAT = 0,
  DOUBLE = 1,
  S64 = 2,
  VEC3_DOUBLE = 3,
  VEC3_FLOAT = 4
};

ControlParameterValueKind getControlParameterValueKindFromString(const std::string& str);

std::string controlParameterValueKindToString(
    ControlParameterValueKind valueKind);

union ControlParameterValuePtr {
  float* f;
  double* d;
  s64* i;
  float* vec3f;
  double* vec3d;
};

union ControlParameterValue {
  float f;
  double d;
  s64 i;
  float vec3f[3];
  double vec3d[3];
};

std::string controlParameterValueToString(ControlParameterValue v,
                                          ControlParameterValueKind kind);

class ControlParameter;

/*!
 * ControlParameterCollections contains a map of all the control parameters.
 */
class ControlParameterCollection {
 public:
  explicit ControlParameterCollection(const std::string& name) : _name(name) {}

  /*!
   * Use this to add a parameter for the first time in the
   * RobotControlParameters or SimulatorControlParameters This should only be
   * used during initialization of a ControlParameter
   */
  void addParameter(ControlParameter* param, const std::string& name) {
    if (mapContains(_map, name)) {
      printf(
          "[ERROR] ControlParameterCollection %s: tried to add parameter %s "
          "twice!\n",
          _name.c_str(), name.c_str());
      throw std::runtime_error("control parameter error");
    }
    _map[name] = param;
  }

  /*!
   * Lookup a control parameter by its name.
   * This does not modify the set field of the control parameter!
   */
  ControlParameter& lookup(const std::string& name) {
    if (mapContains(_map, name)) {
      return *_map[name];
    } else {
      // for now:
      throw std::runtime_error("parameter " + name +
                               " wasn't found in parameter collection " +
                               _name);
    }
  }

  std::string
  printToIniString();  //!< print all control parameters in the INI file format
  std::string printToYamlString();  //!< print all control parameters in the
                                    //!< YAML file format
  bool checkIfAllSet();  //!< are all the control parameters initialized?
  void clearAllSet();
  void deleteAll();

  std::map<std::string, ControlParameter*> _map;

 private:
  std::string _name;
};

class ControlParameter {
 public:
  /*!
   * Constructors for control parameters:
   * Set the type and add to collection.
   * Doesn't "initialize" the parameter - this must be done separately.
   */
  ControlParameter(const std::string& name, double& value,
                   ControlParameterCollection& collection,
                   const std::string& units = "") {
    _name = name;
    truncateName();
    _units = units;
    _value.d = &value;
    _kind = ControlParameterValueKind::DOUBLE;
    collection.addParameter(this, name);
  }

  ControlParameter(const std::string& name, float& value,
                   ControlParameterCollection& collection,
                   const std::string& units = "") {
    _name = name;
    truncateName();
    _units = units;
    _value.f = &value;
    _kind = ControlParameterValueKind::FLOAT;
    collection.addParameter(this, name);
  }

  ControlParameter(const std::string& name, s64& value,
                   ControlParameterCollection& collection,
                   const std::string& units = "") {
    _name = name;
    truncateName();
    _units = units;
    _value.i = &value;
    _kind = ControlParameterValueKind::S64;
    collection.addParameter(this, name);
  }

  ControlParameter(const std::string& name, Vec3<float>& value,
                   ControlParameterCollection& collection,
                   const std::string& units = "") {
    _name = name;
    truncateName();
    _units = units;
    _value.vec3f = value.data();
    _kind = ControlParameterValueKind::VEC3_FLOAT;
    collection.addParameter(this, name);
  }

  ControlParameter(const std::string& name, Vec3<double>& value,
                   ControlParameterCollection& collection,
                   const std::string& units = "") {
    _name = name;
    truncateName();
    _units = units;
    _value.vec3d = value.data();
    _kind = ControlParameterValueKind::VEC3_DOUBLE;
    collection.addParameter(this, name);
  }

  ControlParameter(const std::string& name, ControlParameterValueKind kind) {
    _name = name;
    truncateName();
    _kind = kind;
    _value.vec3d = (double*)&_staticValue;
  }

  void truncateName() {
    if (_name.length() > CONTROL_PARAMETER_MAXIMUM_NAME_LENGTH) {
      printf("[Error] name %s is too long, shortening to ", _name.c_str());
      _name.resize(CONTROL_PARAMETER_MAXIMUM_NAME_LENGTH - 1);
      printf("%s\n", _name.c_str());
    }
  }

  /*!
   * Set initial value of the control parameter.
   * Checks to see that the types are correct
   */
  void initializeDouble(double d) {
    if (_kind != ControlParameterValueKind::DOUBLE) {
      throw std::runtime_error("Tried to initialize control parameter " +
                               _name + " as a double!");
    }
    _set = true;
    *_value.d = d;
  }



  void initializeFloat(float f) {
    if (_kind != ControlParameterValueKind::FLOAT) {
      throw std::runtime_error("Tried to initialize control parameter " +
                               _name + " as a float!");
    }
    _set = true;
    *_value.f = f;
  }

  void initializeInteger(s64 i) {
    if (_kind != ControlParameterValueKind::S64) {
      throw std::runtime_error("Tried to initialize control parameter " +
                               _name + " as an integer!");
    }
    _set = true;
    *_value.i = i;
  }

  void initializeVec3f(const Vec3<float>& v) {
    if (_kind != ControlParameterValueKind::VEC3_FLOAT) {
      throw std::runtime_error("Tried to initialize control parameter " +
                               _name + " as a vector3f");
    }
    _set = true;
    _value.vec3f[0] = v[0];
    _value.vec3f[1] = v[1];
    _value.vec3f[2] = v[2];
  }

  void initializeVec3d(const Vec3<double>& v) {
    if (_kind != ControlParameterValueKind::VEC3_DOUBLE) {
      throw std::runtime_error("Tried to initialize control parameter " +
                               _name + " as a vector3d");
    }
    _set = true;
    _value.vec3d[0] = v[0];
    _value.vec3d[1] = v[1];
    _value.vec3d[2] = v[2];
  }

  void set(ControlParameterValue value, ControlParameterValueKind kind) {
    if (kind != _kind) {
      throw std::runtime_error("type mismatch in set");
    }
    switch (kind) {
      case ControlParameterValueKind::FLOAT:
        *_value.f = value.f;
        break;
      case ControlParameterValueKind::DOUBLE:
        *_value.d = value.d;
        break;
      case ControlParameterValueKind::S64:
        *_value.i = value.i;
        break;
      case ControlParameterValueKind::VEC3_FLOAT:
        _value.vec3f[0] = value.vec3f[0];
        _value.vec3f[1] = value.vec3f[1];
        _value.vec3f[2] = value.vec3f[2];
        break;
      case ControlParameterValueKind::VEC3_DOUBLE:
        _value.vec3d[0] = value.vec3d[0];
        _value.vec3d[1] = value.vec3d[1];
        _value.vec3d[2] = value.vec3d[2];
        break;
      default:
        throw std::runtime_error("invalid kind");
    }
    _set = true;
  }

  ControlParameterValue get(ControlParameterValueKind kind) {
    ControlParameterValue value;
    if (kind != _kind) {
      throw std::runtime_error("type mismatch in get");
    }
    switch (_kind) {
      case ControlParameterValueKind::FLOAT:
        value.f = *_value.f;
        break;
      case ControlParameterValueKind::DOUBLE:
        value.d = *_value.d;
        break;
      case ControlParameterValueKind::S64:
        value.i = *_value.i;
        break;
      case ControlParameterValueKind::VEC3_FLOAT:
        value.vec3f[0] = _value.vec3f[0];
        value.vec3f[1] = _value.vec3f[1];
        value.vec3f[2] = _value.vec3f[2];
        break;
      case ControlParameterValueKind::VEC3_DOUBLE:
        value.vec3d[0] = _value.vec3d[0];
        value.vec3d[1] = _value.vec3d[1];
        value.vec3d[2] = _value.vec3d[2];
        break;
      default:
        throw std::runtime_error("invalid kind");
    }
    return value;
  }

  /*!
   * Convert the value to a string that works in a YAML file
   */
  std::string toString() {
    std::string result;
    switch (_kind) {
      case ControlParameterValueKind::DOUBLE:
        result += numberToString(*_value.d);
        break;
      case ControlParameterValueKind::FLOAT:
        result += numberToString(*_value.f);
        break;
      case ControlParameterValueKind::S64:
        result += std::to_string(*_value.i);
        break;
      case ControlParameterValueKind::VEC3_FLOAT:
        result += "[";
        result += numberToString(_value.vec3f[0]) + ", ";
        result += numberToString(_value.vec3f[1]) + ", ";
        result += numberToString(_value.vec3f[2]) + "]";
        break;
      case ControlParameterValueKind::VEC3_DOUBLE:
        result += "[";
        result += numberToString(_value.vec3d[0]) + ", ";
        result += numberToString(_value.vec3d[1]) + ", ";
        result += numberToString(_value.vec3d[2]) + "]";
        break;
      default:
        result += "<unknown type " + std::to_string((u32)(_kind)) +
                  "> (add it yourself in ControlParameters.h!)";
        break;
    }
    return result;
  }

  bool setFromString(const std::string& value) {
    switch (_kind) {
      case ControlParameterValueKind::DOUBLE:
        *_value.d = std::stod(value);
        break;
      case ControlParameterValueKind::FLOAT:
        *_value.f = std::stof(value);
        break;
      case ControlParameterValueKind::S64:
        *_value.i = std::stoll(value);
        break;
      case ControlParameterValueKind::VEC3_FLOAT: {
        Vec3<float> v = stringToVec3<float>(value);
        _value.vec3f[0] = v[0];
        _value.vec3f[1] = v[1];
        _value.vec3f[2] = v[2];
      } break;
      case ControlParameterValueKind::VEC3_DOUBLE: {
        Vec3<double> v = stringToVec3<double>(value);
        _value.vec3d[0] = v[0];
        _value.vec3d[1] = v[1];
        _value.vec3d[2] = v[2];
      } break;
      default:
        return false;
    }
    return true;
  }

  bool _set = false;
  ControlParameterValuePtr _value;
  std::string _name;
  std::string _units;
  ControlParameterValueKind _kind;
  ControlParameterValue _staticValue;

 private:
};

/*!
 * Parent class for groups of parameters
 * RobotParameters and SimulatorParameters inherit from this class
 */
class ControlParameters {
 public:
  /*!
   * Each control parameter group must have a unique name so the ini files don't
   * mixed up
   * @param name
   */
  ControlParameters(const std::string& name) : collection(name), _name(name) {}

  /*!
   * If true, all parameters have been initialized in one way or another
   */
  bool isFullyInitialized() { return collection.checkIfAllSet(); }

  /*!
   * Directly initialize a given control parameter
   */
  void initializeDouble(const std::string& name, double d) {
    collection.lookup(name).initializeDouble(d);
  }

  void initializeFloat(const std::string& name, float f) {
    collection.lookup(name).initializeFloat(f);
  }

  void initializeInteger(const std::string& name, s64 i) {
    collection.lookup(name).initializeInteger(i);
  }

  void initializeVec3f(const std::string& name, Vec3<float>& v) {
    collection.lookup(name).initializeVec3f(v);
  }

  void initializeVec3d(const std::string& name, Vec3<double>& v) {
    collection.lookup(name).initializeVec3d(v);
  }

  void lockMutex() { _mutex.lock(); }

  void unlockMutex() { _mutex.unlock(); }

  void writeToIniFile(const std::string& path);
  void initializeFromIniFile(const std::string& path);
  void initializeFromYamlFile(const std::string& path);
  void defineAndInitializeFromYamlFile(const std::string& path);

  void writeToYamlFile(const std::string& path);

  std::string generateUnitializedList();

  ControlParameterCollection collection;

 protected:
  std::string _name;
  std::mutex _mutex;
};

#endif  // PROJECT_CONTROLPAREMETERS_H
