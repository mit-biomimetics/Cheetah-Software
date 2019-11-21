#!/bin/bash


GREEN='\033[0;32m'
NC='\033[0m' # No Color

echo -e "${GREEN} Starting LCM type generation...${NC}"

cd ../lcm-types
# Clean
rm */*.jar
rm */*.java
rm */*.hpp
rm */*.class
rm */*.py
rm */*.pyc

# Make
lcm-gen -xp *.lcm
mkdir -p cpp
mv *.hpp cpp

mkdir -p python
mv *.py python

FILES=$(ls */*.class)
echo ${FILES} > file_list.txt


echo -e "${GREEN} Done with LCM type generation${NC}"
