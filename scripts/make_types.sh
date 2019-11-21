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
lcm-gen -jxp *.lcm
cp /usr/local/share/java/lcm.jar .
javac -cp lcm.jar */*.java
jar cf my_types.jar */*.class
mkdir -p java
mv my_types.jar java
mv lcm.jar java
mkdir -p cpp
mv *.hpp cpp

mkdir -p python
mv *.py python

FILES=$(ls */*.class)
echo ${FILES} > file_list.txt


echo -e "${GREEN} Done with LCM type generation${NC}"
