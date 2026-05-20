#ifndef _CONFIG_H_
#define _CONFIG_H_

// 1. Total packing area (L * L)
const double AREA = 10000.0;

// 2. Coefficient for dimensionless packing time (rsa_step = coefficient * area)
const double RSA_STEP_COEFFICIENT = 1E4;

// 3. Minimum number of trials per voxel in final refinement stage (second parameter of final_addnum)
const double FINAL_ADDNUM_MIN = 1000.0;

// 4. Number of independent simulation cycles
const int CYCLE_NUM = 1;

#endif // _CONFIG_H_
