
#cp $CHEETAH_SOFTWARE_PATH/Analysis/PR-MPC/CodeGenPRMPC/PRMPC/codegen/lib/MPCBounds/MPC*.cpp .
#cp $CHEETAH_SOFTWARE_PATH/Analysis/PR-MPC/CodeGenPRMPC/PRMPC/codegen/lib/MPCBounds/MPC*.h .

#rm *_*.cpp
#rm *_*.h

#for file in *.cpp
#do
    #sed -i '/#include "rt_nonfinite.h"/d' "$file"
#done

#for file in *.h
#do
    #sed -i '/#include "rt_nonfinite.h"/d' "$file"
    #sed -i '/#include "rtwtypes.h"/d' "$file"
    #sed -i '/#include "MPCBounds_types.h"/d' "$file"
#done

for file in *SP.*
do
    sed -i -e 's/double/int/g' "$file"
    sed -i -e 's/\.0//g' "$file"
done

#sed -i -e 's/rt_powd_snf/pow/g' ./MPCInitialize.cpp

