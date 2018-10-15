#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <iterator>
#include <ctime>
#include <unordered_map>
#include <algorithm>
#include <cmath>
#include <chrono>
#include <iostream>
#include <fstream>
#include <string>
#include <iterator>
#include <unordered_map>
#include <set>
#include <algorithm>
#include <chrono>
#include <map>
#include <set>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <pthread.h>
#include "bmean.h"
#include "utils.h"



using namespace std;
using namespace chrono;




int main(int argc, char ** argv){
	vector<string> input;



	//~ input.push_back("AAAATCCCCTTTGGGG");
	//~ input.push_back("AAAAGCCCCCCCGGGG");
	//~ input.push_back("AAAAGCACCAAAGGGG");
	//~ auto result(MSABMAAC(input,4));



	string str(rand_seq(1000));
	input.push_back(str);
	for(uint i(0);i<60;++i){
		string str_mut(str);
		mutate(str_mut,str_mut.size()*10/100);
		input.push_back(str_mut);
	}
	cout<<input[5]<<endl;
	auto result(MSABMAAC(input,8));



	//~ for(uint32_t iR(0);iR<result.size();++iR){
		//~ cout<<input[iR]<<" splitted in:			";
		//~ for(uint32_t is(0);is<result[iR].size();++is){
			//~ cout<<result[iR][is]<<"	|";
		//~ }
		//~ cout<<endl;
	//~ }
	//~ cout<<result.size()<<endl<<endl;
	for(uint32_t is(0);is<result[0].size();++is){
		cout<<"-----------------------------------------------------------------------------------------------------------------------------------------"<<endl;
		for(uint32_t iR(0);iR<result.size();++iR){
			cout<<result[iR][is]<<endl;
		}
	}
	return 0;
}
