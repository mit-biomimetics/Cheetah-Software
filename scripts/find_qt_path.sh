#!/bin/bash

# 优先检查/opt/qt5.12
if [ -d "/opt/qt5.12" ]; then
    printf "/opt/qt5.12"
    exit 0
fi

# 兼容原有的~/Qt目录
if [ -d "${HOME}/Qt/" ]; then
    QT_VER="$(ls ~/Qt/ | grep 5 -m1)"
    printf "${HOME}/Qt/${QT_VER}/gcc_64/"
    exit 0
fi

# 如果都找不到，输出空字符串
printf ""