@echo on


:: Perform C Tests with standard configuration
:: -----------------------------------------------------
cd %APPVEYOR_BUILD_FOLDER%
mkdir build
cd build
cmake -G "%CMAKE_PROJECT%" -DUNITTESTS=ON ..
cmake --build .
%APPVEYOR_BUILD_FOLDER%\build\out\osqp_tester.exe
if errorlevel 1 exit /b 1

:: Perform C Tests with floats
:: -----------------------------------------------------
cd %APPVEYOR_BUILD_FOLDER%
rmdir /s /q build
mkdir build
cd build
cmake -G "%CMAKE_PROJECT%" -DDFLOAT=ON -DUNITTESTS=ON ..
cmake --build .
%APPVEYOR_BUILD_FOLDER%\build\out\osqp_tester.exe
if errorlevel 1 exit /b 1

:: Perform C Tests with short integers
:: -----------------------------------------------------
cd %APPVEYOR_BUILD_FOLDER%
rmdir /s /q build
mkdir build
cd build
cmake -G "%CMAKE_PROJECT%" -DDLONG=OFF -DUNITTESTS=ON ..
cmake --build .
%APPVEYOR_BUILD_FOLDER%\build\out\osqp_tester.exe
if errorlevel 1 exit /b 1

@echo off
