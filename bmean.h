#ifndef BMEAN
#define BMEAN



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




using namespace std;





vector<vector<string>> MSABMAAC(const vector<string>& nadine,uint32_t la,double cuisine);



#endif
