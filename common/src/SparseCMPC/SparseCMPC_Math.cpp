#include "SparseCMPC/SparseCMPC.h"
#include <eigen3/Eigen/Dense>
#include <eigen3/unsupported/Eigen/MatrixFunctions>


void SparseCMPC::c2d(u32 trajIdx, u32 bBlockStartIdx, u32 block_count) {

  Eigen::Matrix<double, 24, 24> AB, expmm;
  AB.setZero();
  AB.block(0,0,12,12) = _aMat[trajIdx];

  for(u32 i = bBlockStartIdx; i < bBlockStartIdx + block_count; i++) {
    BblockID id = _bBlockIds[i];
    if(id.timestep != trajIdx) throw std::runtime_error("c2d timestep error");
    AB.block(0,12 + 3 * id.foot, 12, 3) = _bBlocks[i];
  }

  AB *= _dtTrajectory[trajIdx];

  expmm = AB.exp();
  _aMat[trajIdx] = expmm.block(0,0,12,12);

  for(u32 i = bBlockStartIdx; i < bBlockStartIdx + block_count; i++) {
    //BblockID id = _bBlockIds[i];
    //_bBlocks[i] = expmm.block(0,12 + 3 * id.foot, 12, 3);
    _bBlocks[i] *= _dtTrajectory[trajIdx];
  }
}
