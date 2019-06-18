@echo on

REM Needed to enable to define OSQP_DEPLOY_DIR within the file
@setlocal enabledelayedexpansion

IF "%APPVEYOR_REPO_TAG%" == "true" (

    REM Build C libraries
    cd %APPVEYOR_BUILD_FOLDER%
    del /F /Q build
    mkdir build
    cd build
    cmake -G "%CMAKE_PROJECT%" ..
    cmake --build .

    REM Go to output folder
    cd %APPVEYOR_BUILD_FOLDER%\build\out

    IF "%PLATFORM%" == "x86" (
        set "OSQP_DEPLOY_DIR=osqp-0.5.0-windows32"
    ) ELSE (
        set "OSQP_DEPLOY_DIR=osqp-0.5.0-windows64"
    )
    REM Create directories
    REM NB. We force expansion of the variable at execution time!
    mkdir !OSQP_DEPLOY_DIR!
    mkdir !OSQP_DEPLOY_DIR!\lib
    mkdir !OSQP_DEPLOY_DIR!\include

    REM Copy License
    xcopy ..\..\LICENSE !OSQP_DEPLOY_DIR!

    REM Copy includes
    xcopy ..\..\include\*.h !OSQP_DEPLOY_DIR!\include

    REM Copy static library
    xcopy libosqp.a !OSQP_DEPLOY_DIR!\lib

    REM Copy shared library
    xcopy libosqp.dll !OSQP_DEPLOY_DIR!\lib

    REM Compress package
    7z a -ttar !OSQP_DEPLOY_DIR!.tar !OSQP_DEPLOY_DIR!
    7z a -tgzip !OSQP_DEPLOY_DIR!.tar.gz !OSQP_DEPLOY_DIR!.tar

    REM Deploy to Bintray
    curl -T !OSQP_DEPLOY_DIR!.tar.gz -ubstellato:%BINTRAY_API_KEY% -H "X-Bintray-Package:OSQP" -H "X-Bintray-Version:0.5.0" -H "X-Bintray-Override: 1" https://api.bintray.com/content/bstellato/generic/OSQP/0.5.0/
    if errorlevel 1 exit /b 1

    REM Publish
    curl -X POST -ubstellato:%BINTRAY_API_KEY% https://api.bintray.com/content/bstellato/generic/OSQP/0.5.0/publish
    if errorlevel 1 exit /b 1



REM Close parenthesis for deploying only if it is a tagged commit
)
@echo off
