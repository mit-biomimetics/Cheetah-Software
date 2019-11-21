#!/bin/bash
set -ev

OSQP_VERSION="0.5.0"

# Update variables from install
# CMake
if [[ "${TRAVIS_OS_NAME}" == "linux" ]]; then
    export PATH=${DEPS_DIR}/cmake/bin:${PATH}
fi
# Anaconda
export PATH=${DEPS_DIR}/miniconda/bin:$PATH
hash -r
source activate testenv


echo "Creating Bintray package..."

# Compile OSQP
cd ${TRAVIS_BUILD_DIR}
rm -rf build
mkdir build
cd build
cmake -G "Unix Makefiles" ..
make

cd ${TRAVIS_BUILD_DIR}/build/out
if [[ "$TRAVIS_OS_NAME" == "osx" ]]; then
    OS_NAME="mac"
    OS_SHARED_LIB_EXT="dylib"
else
    OS_NAME="linux"
    OS_SHARED_LIB_EXT="so"
fi
OSQP_DEPLOY_DIR=osqp-${OSQP_VERSION}-${OS_NAME}64
mkdir $OSQP_DEPLOY_DIR/
mkdir $OSQP_DEPLOY_DIR/lib
mkdir $OSQP_DEPLOY_DIR/include
# Copy license
cp ../../LICENSE $OSQP_DEPLOY_DIR/
# Copy includes
cp ../../include/*.h  $OSQP_DEPLOY_DIR/include
# Copy static library
cp libosqp.a $OSQP_DEPLOY_DIR/lib
# Copy shared library
cp libosqp.$OS_SHARED_LIB_EXT $OSQP_DEPLOY_DIR/lib
# Compress package
tar -czvf $OSQP_DEPLOY_DIR.tar.gz  $OSQP_DEPLOY_DIR


# Deploy package
curl -T $OSQP_DEPLOY_DIR.tar.gz -ubstellato:$BINTRAY_API_KEY -H "X-Bintray-Package:OSQP" -H "X-Bintray-Version:${OSQP_VERSION}" -H "X-Bintray-Override: 1" https://api.bintray.com/content/bstellato/generic/OSQP/${OSQP_VERSION}/


echo "Creating Bintray sources package..."

OSQP_SOURCES=osqp-${OSQP_VERSION}

# Clone OSQP repository
cd ${TRAVIS_BUILD_DIR}
mkdir sources/
cd sources/
git clone https://github.com/$TRAVIS_REPO_SLUG.git ${OSQP_SOURCES} --recursive
cd ${OSQP_SOURCES}
git checkout -qf $TRAVIS_COMMIT
git submodule update
cd ..

# Create archive ignoring hidden files
tar --exclude=".*" -czvf ${OSQP_SOURCES}.tar.gz ${OSQP_SOURCES}

# Deploy sources
curl -T ${OSQP_SOURCES}.tar.gz -ubstellato:$BINTRAY_API_KEY -H "X-Bintray-Package:OSQP" -H "X-Bintray-Version:${OSQP_VERSION}" -H "X-Bintray-Override: 1" https://api.bintray.com/content/bstellato/generic/OSQP/${OSQP_VERSION}/


# Publish deployed files
curl -X POST -ubstellato:$BINTRAY_API_KEY https://api.bintray.com/content/bstellato/generic/OSQP/0.5.0/publish


exit 0
