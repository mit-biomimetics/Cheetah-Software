#include "ParamHandler.hpp"
#include <iostream>

ParamHandler::ParamHandler(const std::string &file_name) {
  try {
    config_ = dynacore_YAML::LoadFile(file_name);
    fileLoaded = true;
  } catch (std::exception& e) {
    fileLoaded = false;
  }
}

ParamHandler::~ParamHandler() {}

bool ParamHandler::getString(const std::string &key, std::string &str_value) {
  try {
    str_value = config_[key].as<std::string>();
  } catch (std::exception &e) {
    return false;
  }
  return true;
}

bool ParamHandler::getString(const std::string &category, const std::string &key, std::string &str_value) {
  try {
    str_value = config_[category][key].as<std::string>();
  } catch(std::exception &e) {
    return false;
  }
  return true;
}

bool ParamHandler::getBoolean(const std::string &key, bool &bool_value) {
  try {
    bool_value = config_[key].as<bool>();
  } catch (std::exception &e) {
    return false;
  }
  return true;
}

bool ParamHandler::getInteger(const std::string &key, int &int_value) {
  try {
    int_value = config_[key].as<int>();
  } catch (std::exception &e) {
    return false;
  }
  return true;
}
