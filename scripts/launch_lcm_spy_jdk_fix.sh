#!/bin/bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
cd ${DIR}/../lcm-types/java
export CLASSPATH=${DIR}/../lcm-types/java/my_types.jar
pwd
exec java -server -Djava.net.preferIPv4Stack=true -Xmx128m -Xms64m -ea -cp /usr/local/share/java/lcm.jar lcm.spy.Spy
