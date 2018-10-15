#include <stdio.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <atomic>
#include <mutex>
#include <stdint.h>
#include <unordered_map>
#include <algorithm>
#include "utils.h"



using namespace std;



struct localisation {
	uint32_t read_id;
	int32_t position;
};



struct score_chain {
	int32_t length;
	int32_t score;
	int32_t next_anchor;
};



typedef unordered_map<kmer,vector<localisation>> kmer2localisation;



void fill_index_kmers(const vector<string>& Reads,kmer2localisation& kmer_index,uint32_t kmer_size){
	string read;
	uint32_t offsetUpdateKmer=1<<(2*kmer_size);
	unordered_map<kmer,bool> repeated_kmer;
	localisation here({0,0});
	for(uint32_t iR(0);iR<Reads.size();++iR){
		unordered_map<kmer,uint32_t> local_kmer;
		here.position=0;
		here.read_id=iR;
		read=Reads[iR];
		if(read.size()<kmer_size){continue;}
		kmer seq(str2num(read.substr(0,kmer_size)));
		kmer_index[seq].push_back(here);
		if(++local_kmer[seq]>1){
			repeated_kmer[seq]=true;
		}
		for(uint32_t ir(0);kmer_size+ir<read.size();++ir){
			updateK(seq,read[kmer_size+ir],offsetUpdateKmer);
			++here.position;
			kmer_index[seq].push_back(here);
			if(++local_kmer[seq]>1){
				repeated_kmer[seq]=true;
			}
		}
	}
	auto it = repeated_kmer.begin();
	while(it != repeated_kmer.end()){
		kmer_index.erase(it->first);
		++it;
	}
}



unordered_map<kmer,uint32_t> filter_index_kmers(kmer2localisation& kmer_index, double amount){
	unordered_map<kmer,uint32_t> result;
	//~ cout<<"kmer INdex size before cleaning"<<kmer_index.size()<<endl;
	vector<kmer> to_suppress;
	vector<uint32_t> read_ids;
	auto it = kmer_index.begin();
	while(it != kmer_index.end()){
		for(uint32_t i(0);i<it->second.size();++i){
			read_ids.push_back(it->second[i].read_id);
		}
		//AVOID TO COUNT MULTIPLE OCCURENCE OF A KMER WITHIN A READ
		sort( read_ids.begin(), read_ids.end() );
		int uniqueCount = unique(read_ids.begin(), read_ids.end()) - read_ids.begin();
		if(uniqueCount<amount){
			to_suppress.push_back(it->first);
		}else{
			result[it->first]=uniqueCount;
		}
		++it;
		read_ids.clear();
	}
	for(uint32_t i(0);i<to_suppress.size();++i){
		kmer_index.erase(to_suppress[i]);
	}
	//~ cout<<"kmer INdex size after cleaning"<<kmer_index.size()<<endl;
	return result;
}



bool order_according2read_id (localisation i,localisation j) { return (i.read_id<j.read_id); }



int anchors_ordered_according2reads(const kmer kmer1,const kmer kmer2,  kmer2localisation& kmer_index){
	uint32_t result(0);
	auto v_loc1(kmer_index[kmer1]);
	auto v_loc2(kmer_index[kmer2]);
	//~ sort (v_loc1.begin(), v_loc1.end(), order_according2read_id);
	//~ sort (v_loc2.begin(), v_loc2.end(), order_according2read_id);
	uint32_t i1(0),i2(0);
	//BIG QUESTION HOW TO HANDLE REPEATED KMER HERE
	while(i1<v_loc1.size() and i2<v_loc2.size()){
		if(v_loc1[i1].read_id==v_loc2[i2].read_id){
			if(v_loc1[i1].position>v_loc2[i2].position){
				return -1;
				//COULD ADD A NO IF POSITIONS ARE TOO FAR LIKE IN MINIMAP
			}else{
				++i1;
				++i2;
				++result;
			}
		}else if(v_loc1[i1].read_id<v_loc2[i2].read_id){
			i1++;
		}else{
			i2++;
		}
	}
	return result;
}



score_chain longest_ordered_chain_from_anchors( kmer2localisation& kmer_index, unordered_map<uint,score_chain>& best_chain_computed, uint32_t start, const vector<kmer>& template_read){
	if(best_chain_computed.count(start)==1){
		return best_chain_computed[start];
	}
	int32_t max_chain(-1),max_score(0);
	int32_t next_anchor(-1);
	for(uint i(start+1);i<template_read.size();++i){
		kmer next(template_read[i]);
		int score(anchors_ordered_according2reads(template_read[start],next,kmer_index));
		if(score>0){

			auto p=longest_ordered_chain_from_anchors(kmer_index,best_chain_computed,i,template_read);
			if(p.length>max_chain){
				max_chain=p.length;
				max_score=p.score+score;
				next_anchor=i;
			}else if(p.length==max_chain and p.score+score>max_score) {
				max_score=p.score+score;
				next_anchor=i;
			}else{
			}
		}
	}

	//~ cout<<"SCORE of "<<start<<": "<<max_chain+1<<" "<<max_score<<" "<<next_anchor<<endl;
	best_chain_computed[start]={max_chain+1,max_score,next_anchor};
	return {max_chain+1,max_score,next_anchor};
}



vector<kmer> get_template( kmer2localisation& kmer_index,const string& read,int kmer_size){
	vector<kmer> result;
	uint32_t offsetUpdateKmer=1<<(2*kmer_size);
	kmer seq(str2num(read.substr(0,kmer_size)));
	if(kmer_index.count(seq)){
		result.push_back(seq);
	}
	for(uint32_t ir(0);kmer_size+ir<read.size();++ir){
		updateK(seq,read[kmer_size+ir],offsetUpdateKmer);
		if(kmer_index.count(seq)){
			result.push_back(seq);
		}
	}
	return result;
}




vector<kmer> longest_ordered_chain( kmer2localisation& kmer_index,const vector<kmer>& template_read){
	unordered_map<uint,score_chain> best_chain_computed;
    vector<kmer> result;
    int32_t max_chain(0),max_score(0);
	int32_t next_anchor(-1);
    for(int32_t i(template_read.size()-1);i>=0;--i){
        auto p=longest_ordered_chain_from_anchors(kmer_index,best_chain_computed,i,template_read);
        if(p.length>max_chain){
			max_chain=p.length;
			max_score=p.score;
			next_anchor=i;
		}else if(p.length==max_chain and p.score>max_score) {
			max_score=p.score;
			next_anchor=i;
		}
    }
    while(next_anchor!=-1){
        result.push_back(template_read[next_anchor]);
        next_anchor=best_chain_computed[next_anchor].next_anchor;
    }
	return result;
}



vector<double> average_distance_next_anchor(kmer2localisation& kmer_index,  vector<kmer>& anchors,unordered_map<kmer,uint32_t>& k_count, bool clean){
	vector<double> result;
	vector<uint32_t> v_sum;
	vector<uint32_t> v_count;
	vector<kmer> curated_anchors;
	uint32_t min_distance(5);

	for(uint i(0);i+1<anchors.size();++i){
		uint32_t sum(0),count(0);
		auto v_loc1(kmer_index[anchors[i]]);//THEY SHOULD BE READS SORTED
		auto v_loc2(kmer_index[anchors[i+1]]);
		//~ sort (v_loc1.begin(), v_loc1.end(), order_according2read_id);
		//~ sort (v_loc2.begin(), v_loc2.end(), order_according2read_id);
		uint32_t i1(0),i2(0);
		while(i1<v_loc1.size() and i2<v_loc2.size()){
			if(v_loc1[i1].read_id==v_loc2[i2].read_id){
				sum+=v_loc2[i2].position-v_loc1[i1].position;
				++count;
				++i1;
				++i2;
			}else if(v_loc1[i1].read_id<v_loc2[i2].read_id){
				i1++;
			}else{
				i2++;
			}
		}
		if(count!=0){
			v_sum.push_back(sum);
			v_count.push_back(count);

			result.push_back(sum/count);
		}else{
			cout<<"SHOULD NOT HAPPEND"<<endl;cin.get();
			result.push_back(-1);
		}
	}

	if(clean){
		for(uint32_t i(0);i+1<anchors.size();++i){
			if(result[i]<min_distance){
				kmer heaviest_anchor(anchors[i]);
				uint32_t weight(k_count[anchors[i]]);
				double new_distance(0);
				while(i+1<anchors.size() and new_distance<min_distance){
					++i;
					if(weight<=k_count[anchors[i]]){
						weight=k_count[anchors[i]];
						heaviest_anchor=(anchors[i]);
					}
					new_distance+=result[i-1];
				}
				curated_anchors.push_back(heaviest_anchor);

			}else{
				curated_anchors.push_back(anchors[i]);
			}
		}
		anchors=curated_anchors;
	}


	return result;
}



int32_t get_position(kmer2localisation& kmer_index,kmer query, uint32_t read_id){
	auto V(kmer_index[query]);
	for(uint32_t i(0);i<V.size();++i){
		if(V[i].read_id==read_id){
			return V[i].position;
		}
	}
	return -1;
}



vector<vector<string>> split_reads(const vector<kmer>& anchors, const vector<double>& relative_positions, const vector<string>& Reads,  kmer2localisation& kmer_index,uint32_t kmer_size){
	vector<vector<string>> result;
	for(uint32_t iR(0);iR<Reads.size();++iR){
		string read=Reads[iR];
		vector<string> split(anchors.size()+1);
		//FIRST AND LAST REGION
		int32_t anchor_position(get_position(kmer_index,anchors[0],iR));
		if(anchor_position!=-1){
			split[0]=read.substr(0,anchor_position);
		}
		anchor_position=(get_position(kmer_index,anchors[anchors.size()-1],iR));
		if(anchor_position!=-1){
			split[anchors.size()]=read.substr(anchor_position);
		}

		for(uint32_t iA(0);iA+1<anchors.size();++iA){
			int32_t anchor_position1(get_position(kmer_index,anchors[iA],iR));
			int32_t anchor_position2(get_position(kmer_index,anchors[iA+1],iR));
			if(anchor_position1!=-1){
				if(anchor_position2!=-1){
					//REGION WITH BOtH ANCHORS
					split[iA+1]=read.substr(anchor_position1,anchor_position2-anchor_position1);
				}else{
					//GOT THE LEFT ANCHOR
					split[iA+1]=read.substr(anchor_position1,relative_positions[iA]);
				}
			}else{
				if(anchor_position2!=-1){
					//GOT THE RIGHT ANCHOR
					if(anchor_position2>relative_positions[iA]){
						split[iA+1]=read.substr(anchor_position2-relative_positions[iA],relative_positions[iA]);
					}
				}
			}
		}
		result.push_back(split);
	}
	return result;
}



vector<vector<string>> MSABMAAC(const vector<string>& Reads,uint32_t k){
	int kmer_size(k);
	double percent_shared(0.6);


	kmer2localisation kmer_index;
	fill_index_kmers(Reads,kmer_index,kmer_size);
	//~ cout<<"PHASE 1 done"<<endl;

	auto kmer_count(filter_index_kmers(kmer_index,percent_shared*(double)Reads.size()));
	//~ cout<<"PHASE 2.1 done"<<endl;
	auto template_read(get_template(kmer_index,Reads[0],kmer_size));
	//~ cout<<"PHASE 2 done"<<endl;

	vector<kmer> anchors(longest_ordered_chain(kmer_index, template_read));
	//~ cout<<"PHASE 3 done"<<endl;

	//~ vector<double> relative_positions(average_distance_next_anchor(kmer_index,anchors,kmer_count,true));
	//~ relative_positions=(average_distance_next_anchor(kmer_index,anchors,kmer_count,true));
	vector<double> relative_positions=(average_distance_next_anchor(kmer_index,anchors,kmer_count,false));
	//~ cout<<"PHASE 4 done"<<endl;

	vector<vector<string>> result(split_reads(anchors,relative_positions,Reads,kmer_index,kmer_size));
	//~ cout<<"PHASE 5 done"<<endl;

	return result;
}
