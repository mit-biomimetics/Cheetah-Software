#!/bin/bash
set -ev

# Create deps dir
mkdir ${DEPS_DIR}
cd ${DEPS_DIR}

# Install CMake
if [[ "${TRAVIS_OS_NAME}" == "linux" ]]; then
    CMAKE_URL="http://www.cmake.org/files/v3.7/cmake-3.7.1-Linux-x86_64.tar.gz"
    mkdir cmake && wget --quiet -O - ${CMAKE_URL} | tar --strip-components=1 -xz -C cmake
    export PATH=${DEPS_DIR}/cmake/bin:${PATH}
else
brew update
brew upgrade cmake || brew install cmake
fi
cmake --version


# Install lcov
if [[ "$TRAVIS_OS_NAME" == "linux" ]]; then
    sudo apt-get install -y lcov
else if [[ "$TRAVIS_OS_NAME" == "osx" ]]; then
    brew update
    brew install lcov
    fi
fi

gem install coveralls-lcov


# Install Anaconda

# Use the miniconda installer for faster download / install of conda
# itself
if [[ $TRAVIS_OS_NAME == 'osx' ]]; then
    wget http://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh -O miniconda.sh;
else
    wget http://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh;
fi
chmod +x miniconda.sh && ./miniconda.sh -b -p ${DEPS_DIR}/miniconda
export PATH=${DEPS_DIR}/miniconda/bin:$PATH
hash -r
conda config --set always_yes yes --set changeps1 no
conda update --yes -q conda
conda create -n testenv --yes python=$PYTHON_VERSION numpy scipy future
source activate testenv
