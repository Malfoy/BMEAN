#ifndef GLOBAL
#define GLOBAL


#include <vector>
#include <string>



extern "C"{
#include "lpo.h"
#include "msa_format.h"
#include "align_score.h"
#include "default.h"
#include "poa.h"
#include "seq_util.h"
}
ResidueScoreMatrix_T score_matrix;
bool score_matrix_init=false;


#endif
