@echo on
:: Make sure all the submodules are updated correctly
cd %APPVEYOR_BUILD_FOLDER%
git submodule update --init --recursive

:: Remove entry with sh.exe from PATH to fix error with MinGW toolchain
:: (For MinGW make to work correctly sh.exe must NOT be in your path)
:: http://stackoverflow.com/a/3870338/2288008
set PATH=%PATH:C:\Program Files\Git\usr\bin;=%


IF "%PLATFORM%"=="x86" (
    set MINGW_PATH=C:\mingw-w64\i686-6.3.0-posix-dwarf-rt_v5-rev1\mingw32\bin
) ELSE (
    set MINGW_PATH=C:\mingw-w64\x86_64-6.3.0-posix-seh-rt_v5-rev1\mingw64\bin
)
set PATH=%MINGW_PATH%;%PATH%


:: Activate test environment anaconda
IF "%PLATFORM%"=="x86" (
	set MINICONDA_PATH=%MINICONDA%
) ELSE (
	set MINICONDA_PATH=%MINICONDA%-%PLATFORM%
)
set PATH=%MINICONDA_PATH%;%MINICONDA_PATH%\\Scripts;%PATH%


conda config --set always_yes yes --set changeps1 no
REM This, together with next line, disables conda auto update (fixes problem with tqdm)
conda config --set auto_update_conda false
REM conda update -q conda
conda info -a
conda create -q -n test-environment python=%PYTHON_VERSION% numpy scipy future
if errorlevel 1 exit /b 1
:: NB: Need to run with call otherwise the script hangs
call activate test-environment


:: Set environment for build if 64bit
:: NB: Needed during conda build!
IF "%PLATFORM%"=="x64" (
call "C:\Program Files\Microsoft SDKs\Windows\v7.1\Bin\SetEnv.cmd" /x64
) ELSE (
REM Set environment for 32bit
call "C:\Program Files (x86)\Microsoft Visual Studio 14.0\VC\vcvarsall.bat" x86
)


@echo off
