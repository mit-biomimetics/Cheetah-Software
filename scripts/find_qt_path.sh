#/bin/bash

QT_VER="$(ls ~/Qt/ | grep 5 -m1)"

printf "${HOME}/Qt/${QT_VER}/gcc_64/"

