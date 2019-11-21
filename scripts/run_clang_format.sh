#!/bin/bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

find ${DIR}/../robot ${DIR}/../sim ${DIR}/../common -name *.h -o -iname *.hpp -o -iname *.cpp -o -iname *.c | xargs clang-format -style=google -i
