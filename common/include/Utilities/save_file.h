#ifndef SAVE_FILE_H
#define SAVE_FILE_H

#include <Configuration.h>
#include <cppTypes.h>
#include <stdio.h>
#include <sys/stat.h>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <list>
#include <string>

static std::list<std::string> gs_fileName_string;  // global & static

void cleaning_file(const std::string& folder_name, const std::string& file_name,
                   std::string& full_ret_file);

void create_folder(const std::string& folder_name);

template <typename T>
void saveVector(const DVec<T>& _vec, const std::string& folder_name,
                const std::string& file_name) {
  std::string full_file_name;
  cleaning_file(folder_name, file_name, full_file_name);

  std::ofstream savefile(full_file_name.c_str(), std::ios::app);
  for (int i(0); i < _vec.rows(); ++i) {
    savefile << _vec(i) << "\t";
  }
  savefile << "\n";
  savefile.flush();
}

template <typename T>
void saveVector(const Vec3<T>& _vec, const std::string& folder_name,
                const std::string& file_name) {
  saveVector((DVec<T>)_vec, folder_name, file_name);
}

template <typename T>
void saveVector(const std::vector<T>& _vec, const std::string& folder_name,
                const std::string& file_name) {
  std::string full_file_name;
  cleaning_file(folder_name, file_name, full_file_name);
  std::ofstream savefile(full_file_name.c_str(), std::ios::app);

  for (unsigned int i(0); i < _vec.size(); ++i) {
    savefile << _vec[i] << "\t";
  }
  savefile << "\n";
  savefile.flush();
}

template <typename T>
void saveVector(T* _vec, const std::string& folder_name,
                const std::string& file_name, int size) {
  std::string full_file_name;
  cleaning_file(folder_name, file_name, full_file_name);
  std::ofstream savefile(full_file_name.c_str(), std::ios::app);
  for (int i(0); i < size; ++i) {
    savefile << _vec[i] << "\t";
  }
  savefile << "\n";
  savefile.flush();
}

template <typename T>
void saveValue(T _value, const std::string& folder_name,
               const std::string& file_name) {
  std::string full_file_name;
  cleaning_file(folder_name, file_name, full_file_name);
  std::ofstream savefile(full_file_name.c_str(), std::ios::app);

  savefile << _value << "\n";
  savefile.flush();
}

#endif
