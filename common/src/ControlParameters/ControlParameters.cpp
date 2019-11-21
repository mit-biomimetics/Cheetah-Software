/*! @file ControlParameters.cpp
 *  @brief Interface to set gains/control parameters for simulator and robot
 *  These are designed to be updated infrequently.  For high frequency data,
 * consider using Driver Inputs or adding to the Robot Debug Data instead.
 */

#include <utility>

#include "ControlParameters/ControlParameters.h"
#include "INIReader.h"
#include "ParamHandler.hpp"
#include "Utilities/utilities.h"


#define YAML_COLLECTION_NAME_KEY "__collection-name__"

/*!
 * Convert control parameter data types into strings
 */
std::string controlParameterValueKindToString(
    ControlParameterValueKind valueKind) {
  switch (valueKind) {
    case ControlParameterValueKind::S64:
      return "s64";
    case ControlParameterValueKind::DOUBLE:
      return "double";
    case ControlParameterValueKind::FLOAT:
      return "float";
    case ControlParameterValueKind::VEC3_DOUBLE:
      return "vec3d";
    case ControlParameterValueKind::VEC3_FLOAT:
      return "vec3f";
    default:
      return "unknown-ControlParameterValueKind";
  }
}

/*!
 * Convert a control parameter value into a human readable string.
 */
std::string controlParameterValueToString(ControlParameterValue value,
                                          ControlParameterValueKind kind) {
  std::string result;
  switch (kind) {
    case ControlParameterValueKind::DOUBLE:
      result += numberToString(value.d);
      break;
    case ControlParameterValueKind::FLOAT:
      result += numberToString(value.f);
      break;
    case ControlParameterValueKind::S64:
      result += std::to_string(value.i);
      break;
    case ControlParameterValueKind::VEC3_FLOAT:
      result += "[";
      result += numberToString(value.vec3f[0]) + ", ";
      result += numberToString(value.vec3f[1]) + ", ";
      result += numberToString(value.vec3f[2]) + "]";
      break;
    case ControlParameterValueKind::VEC3_DOUBLE:
      result += "[";
      result += numberToString(value.vec3d[0]) + ", ";
      result += numberToString(value.vec3d[1]) + ", ";
      result += numberToString(value.vec3d[2]) + "]";
      break;
    default:
      result += "<unknown type " + std::to_string((u32)(kind)) +
                "> (add it yourself in ControlParameterInterface.h!)";
      break;
  }
  return result;
}

/*!
 * Check if all parameters have had their value set
 * @return if all set
 */
bool ControlParameterCollection::checkIfAllSet() {
  for (auto& kv : _map) {
    if (!kv.second->_set) {
      return false;
    }
  }
  return true;
}

/*!
 * Mark all parameters as not set
 */
void ControlParameterCollection::clearAllSet() {
  for (auto& kv : _map) {
    kv.second->_set = false;
  }
}

/*!
 * Remove all parameters
 */
void ControlParameterCollection::deleteAll() {
  for(auto& kv : _map) {
    delete kv.second;
  }
  _map.clear();
}

/*!
 * Write a list of uninitialized parameters to a string
 */
std::string ControlParameters::generateUnitializedList() {
  std::string result;
  for (auto& kv : collection._map) {
    if (!kv.second->_set) {
      result += kv.second->_name + "    :\n";
    }
  }

  return result;
}

/*!
 * Print all parameters to a string that can be used for an INI file
 */
std::string ControlParameterCollection::printToIniString() {
  std::string result = ";; Generated on " + getCurrentTimeAndDate() + "\n";

  // ini section
  result += "[";
  result += _name + "]\n\n";

  std::vector<std::string> lines;

  // names (and max name length)
  int maxLength_name = 0;
  for (auto& kv : _map) {
    maxLength_name = std::max(maxLength_name, (int)kv.first.length());
    lines.push_back(kv.first);
  }

  // pad, equals sign, and number
  size_t i = 0;
  int maxLength_number = 0;
  for (auto& kv : _map) {
    int charsToAdd = maxLength_name - (int)lines[i].length();
    assert(charsToAdd >= 0);
    for (int j = 0; j < charsToAdd; j++) {
      lines[i].push_back(' ');
    }
    lines[i] += " = ";
    lines[i] += kv.second->toString();
    maxLength_number = std::max(maxLength_number, (int)lines[i].length());
    i++;
  }

  // pad, ;;, and units
  i = 0;
  for (auto& kv : _map) {
    int charsToAdd = maxLength_number - (int)lines[i].length();
    assert(charsToAdd >= 0);
    for (int j = 0; j < charsToAdd; j++) {
      lines[i].push_back(' ');
    }
    lines[i] += " ;; ";
    if (kv.second->_units.empty()) {
      lines[i] += "No units specified.  Add them!";
    } else {
      lines[i] += kv.second->_units;
    }

    i++;
  }

  // combine lines
  for (auto& line : lines) {
    result += line + "\n";
  }

  return result;
}

/*!
 * Print all parameters to a string which can be used as a YAML file
 */
std::string ControlParameterCollection::printToYamlString() {
  std::string result = "# Generated on " + getCurrentTimeAndDate() + "\n";

  result += YAML_COLLECTION_NAME_KEY;
  result += ": ";
  result += _name + "\n\n";

  std::vector<std::string> lines;

  // names
  int maxLength_name = 0;
  for (auto& kv : _map) {
    maxLength_name = std::max(maxLength_name, (int)kv.first.length());
    lines.push_back(kv.first);
  }

  // name pad, :, and number
  size_t i = 0;
  int maxLength_number = 0;
  for (auto& kv : _map) {
    int charsToAdd = maxLength_name - (int)lines[i].length();
    assert(charsToAdd >= 0);
    for (int j = 0; j < charsToAdd; j++) {
      lines[i].push_back(' ');
    }
    lines[i] += ": ";
    lines[i] += kv.second->toString();
    maxLength_number = std::max(maxLength_number, (int)lines[i].length());
    i++;
  }

  // combine lines
  for (auto& line : lines) {
    result += line + "\n";
  }

  return result;
}

/*!
 * Write all parameters and values to an INI file
 * @param path : the file name
 */
void ControlParameters::writeToIniFile(const std::string& path) {
  writeStringToFile(path, collection.printToIniString());
}

/*!
 * Write all parameters and values to a YAML file
 * @param path : the file name
 */
void ControlParameters::writeToYamlFile(const std::string& path) {
  writeStringToFile(path, collection.printToYamlString());
}

/*!
 * Load values from an INI file
 * @param path : the file name
 */
void ControlParameters::initializeFromIniFile(const std::string& path) {
  INIReader iniReader(path);
  if (iniReader.ParseError() < 0) {
    printf(
        "[ERROR] Could not open ini file %s : not initializing control "
        "parameters!\n",
        path.c_str());
    throw std::runtime_error("ini file bad");
  }

  std::set<std::string> sections = iniReader.GetSections();

  if (sections.size() != 1) {
    printf(
        "[ERROR] INI file %s had %ld sections (expected 1) : not initializing "
        "control parameters\n",
        path.c_str(), sections.size());
    throw std::runtime_error("ini file bad");
  }

  std::string sectionName = *(sections.begin());

  if (sectionName != _name) {
    printf(
        "[ERROR] INI file %s has section name %s, which cannot be used to "
        "initialize %s\n",
        path.c_str(), sectionName.c_str(), _name.c_str());
    throw std::runtime_error("ini file bad");
  }

  std::set<std::string> parameterNames = iniReader.GetFields(sectionName);

  for (auto& name : parameterNames) {
    ControlParameter& cp = collection.lookup(name);
    switch (cp._kind) {
      case ControlParameterValueKind::DOUBLE:
        cp.initializeDouble(iniReader.GetReal(sectionName, name, 0.));
        break;
      case ControlParameterValueKind::FLOAT:
        cp.initializeFloat((float)iniReader.GetReal(sectionName, name, 0.));
        break;
      case ControlParameterValueKind::S64:
        cp.initializeInteger(iniReader.GetInteger(sectionName, name, 0));
        break;
      default:
        throw std::runtime_error("can't read type " +
                                 std::to_string((u32)cp._kind) +
                                 " from ini file");
        break;
    }
  }
}

/*!
 * Determine the kind of a control parameter from a string representation
 * @param str : the string (like "2", or "[2,3,3]")
 * @return the kind
 */
ControlParameterValueKind getControlParameterValueKindFromString(const std::string& str) {
  if(str.find('[') != std::string::npos) {
    // vec type
    if(str.find('f') != std::string::npos) {
      return ControlParameterValueKind::VEC3_FLOAT;
    } else {
      return ControlParameterValueKind::VEC3_DOUBLE;
    }
  } else {
    // scalar type
    if(str.find('.') != std::string::npos) {
      if(str.find('f') != std::string::npos) {
        return ControlParameterValueKind::FLOAT;
      } else {
        return ControlParameterValueKind::DOUBLE;
      }
    } else {
      // integer
      return ControlParameterValueKind::S64;
    }
  }
}

/*!
 * Add and initialize parameters from a YAML file
 */
void ControlParameters::defineAndInitializeFromYamlFile(const std::string &path) {
  ParamHandler paramHandler(path);

  if (!paramHandler.fileOpenedSuccessfully()) {
    printf(
        "[ERROR] Could not open yaml file %s : not initializing control "
        "parameters!\n",
        path.c_str());
    throw std::runtime_error("yaml file bad");
  }

  std::string name;
  if (!paramHandler.getString(YAML_COLLECTION_NAME_KEY, name)) {
    printf("[ERROR] YAML doesn't have a a collection name field named %s\n",
           YAML_COLLECTION_NAME_KEY);
    throw std::runtime_error("yaml file bad");
  }

  if (name != _name) {
    printf(
        "[ERROR] YAML file %s has collection name %s which cannot be used to "
        "initialize %s\n",
        path.c_str(), name.c_str(), _name.c_str());
    throw std::runtime_error("yaml file bad");
  }

  std::vector<std::string> keys = paramHandler.getKeys();

  for (auto& key : keys) {
    if (key == YAML_COLLECTION_NAME_KEY) continue;
    std::string valueString;
    paramHandler.getString(key, valueString);
    ControlParameterValueKind kind;
    if(valueString.empty()) {
      kind = ControlParameterValueKind::VEC3_DOUBLE;
    } else {
      kind = ControlParameterValueKind::DOUBLE;
    }

    ControlParameter* cp = new ControlParameter(key, kind);
    collection.addParameter(cp, key);
    switch (cp->_kind) {
      case ControlParameterValueKind::DOUBLE: {
        double d;
        assert(paramHandler.getValue(key, d));
        cp->initializeDouble(d);
      } break;

      case ControlParameterValueKind::FLOAT: {
        float f;
        assert(paramHandler.getValue(key, f));
        cp->initializeFloat(f);
      } break;

      case ControlParameterValueKind::S64: {
        s64 f;
        assert(paramHandler.getValue(key, f));
        cp->initializeInteger(f);
      } break;

      case ControlParameterValueKind::VEC3_DOUBLE: {
        std::vector<double> vv;
        assert(paramHandler.getVector(key, vv));
        assert(vv.size() == 3);
        Vec3<double> v(vv[0], vv[1], vv[2]);
        cp->initializeVec3d(v);
      } break;

      case ControlParameterValueKind::VEC3_FLOAT: {
        std::vector<float> vv;
        assert(paramHandler.getVector(key, vv));
        assert(vv.size() == 3);
        Vec3<float> v(vv[0], vv[1], vv[2]);
        cp->initializeVec3f(v);
      } break;

      default:
        throw std::runtime_error("can't read type " +
                                 std::to_string((u32)cp->_kind) +
                                 " from yaml file");
        break;
    }
  }
}

/*!
 * Set parameters from a YAML file
 * @param path : the file name
 */
void ControlParameters::initializeFromYamlFile(const std::string& path) {
  ParamHandler paramHandler(path);

  if (!paramHandler.fileOpenedSuccessfully()) {
    printf(
        "[ERROR] Could not open yaml file %s : not initializing control "
        "parameters!\n",
        path.c_str());
    throw std::runtime_error("yaml file bad");
  }

  std::string name;
  if (!paramHandler.getString(YAML_COLLECTION_NAME_KEY, name)) {
    printf("[ERROR] YAML doesn't have a a collection name field named %s\n",
           YAML_COLLECTION_NAME_KEY);
    throw std::runtime_error("yaml file bad");
  }

  if (name != _name) {
    printf(
        "[ERROR] YAML file %s has collection name %s which cannot be used to "
        "initialize %s\n",
        path.c_str(), name.c_str(), _name.c_str());
    throw std::runtime_error("yaml file bad");
  }

  std::vector<std::string> keys = paramHandler.getKeys();

  for (auto& key : keys) {
    if (key == YAML_COLLECTION_NAME_KEY) continue;
    ControlParameter& cp = collection.lookup(key);
    switch (cp._kind) {
      case ControlParameterValueKind::DOUBLE: {
        double d;
        assert(paramHandler.getValue(key, d));
        cp.initializeDouble(d);
      } break;

      case ControlParameterValueKind::FLOAT: {
        float f;
        assert(paramHandler.getValue(key, f));
        cp.initializeFloat(f);
      } break;

      case ControlParameterValueKind::S64: {
        s64 f;
        assert(paramHandler.getValue(key, f));
        cp.initializeInteger(f);
      } break;

      case ControlParameterValueKind::VEC3_DOUBLE: {
        std::vector<double> vv;
        assert(paramHandler.getVector(key, vv));
        assert(vv.size() == 3);
        Vec3<double> v(vv[0], vv[1], vv[2]);
        cp.initializeVec3d(v);
      } break;

      case ControlParameterValueKind::VEC3_FLOAT: {
        std::vector<float> vv;
        assert(paramHandler.getVector(key, vv));
        assert(vv.size() == 3);
        Vec3<float> v(vv[0], vv[1], vv[2]);
        cp.initializeVec3f(v);
      } break;

      default:
        throw std::runtime_error("can't read type " +
                                 std::to_string((u32)cp._kind) +
                                 " from yaml file");
        break;
    }
  }
}
