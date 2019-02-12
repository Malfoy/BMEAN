#ifndef BMEAN
#define BMEAN



#include <vector>
#include <string>
#include "utils.h"
#include <cmath>


extern "C"{
#include "lpo.h"
#include "msa_format.h"
#include "align_score.h"
#include "default.h"
#include "poa.h"
#include "seq_util.h"
}




using namespace std;





std::pair<std::vector<std::vector<std::string>>, std::unordered_map<kmer, unsigned>> MSABMAAC(const vector<string>& nadine,uint32_t la,double cuisine, unsigned solidThresh, unsigned minAnchors, unsigned maxMSA);



#endif
