#!/bin/bash

all_num=(2 3 4 5 6 7 8 9 10 11 12 14 15 16 17 18 19 20)

for groupnum in ${all_num[@]}; do
cp "./KPP_Fun/main_"$groupnum".F90" main.F90
./compile
done
exit 0
#EOC
