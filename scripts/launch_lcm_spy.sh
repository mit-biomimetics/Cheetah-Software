#!/bin/bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
cd ${DIR}/../lcm-types/java
export CLASSPATH=${DIR}/../lcm-types/java/my_types.jar
pwd
lcm-spy
