#include <stdio.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <set>
#include <atomic>
#include <mutex>
#include <stdint.h>
#include <unordered_map>
#include <algorithm>
#include <math.h>
#include "utils.h"
#include "bmean.h"
#include "Complete-Striped-Smith-Waterman-Library/src/ssw_cpp.h"
#include "global.h"

extern "C"{
#include "lpo.h"
#include "msa_format.h"
#include "align_score.h"
#include "default.h"
#include "poa.h"
#include "seq_util.h"
}



//TODO REMOVE HASHTABLES

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



void fill_index_kmers(const vector<string>& Reads,kmer2localisation& kmer_index,uint32_t kmer_size, std::unordered_map<kmer, unsigned>& merCounts, unsigned solidThresh){
	std::unordered_map<kmer, unsigned> tmpMerCounts;
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
		// tmpMerCounts[seq]++;
		if(++local_kmer[seq]>1){
			repeated_kmer[seq]=true;
		}
		for(uint32_t ir(0);kmer_size+ir<read.size();++ir){
			updateK(seq,read[kmer_size+ir],offsetUpdateKmer);
			++here.position;
			kmer_index[seq].push_back(here);
			// tmpMerCounts[seq]++;
			if(++local_kmer[seq]>1){
				repeated_kmer[seq]=true;
			}
		}
	}
	// auto it = repeated_kmer.begin();
	// while(it != repeated_kmer.end()){
	// 	kmer_index.erase(it->first);
	// 	++it;
	// }

	// for (auto p : tmpMerCounts) {
	// 	if (p.second >= solidThresh) {
	// 		merCounts[p.first] = p.second;
	// 	}
	// }

	for (auto p : kmer_index) {
		if (p.second.size() >= solidThresh) {
			merCounts[p.first] = p.second.size();
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
	//~ cerr<<"kmer INdex size before cleaning"<<kmer_index.size()<<endl;
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
	//~ cerr<<"kmer INdex size after cleaning"<<kmer_index.size()<<endl;
	return result;
}



bool order_according2read_id (localisation i,localisation j) { return (i.read_id<j.read_id); }



int anchors_ordered_according2reads(const kmer kmer1,const kmer kmer2,  kmer2localisation& kmer_index){
	int32_t result(0);
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



score_chain longest_ordered_chain_from_anchors( kmer2localisation& kmer_index, unordered_map<uint,score_chain>& best_chain_computed, uint32_t start, const vector<kmer>& template_read,double edge_solidity){
	if(best_chain_computed.count(start)==1){
		return best_chain_computed[start];
	}
	int32_t max_chain(-1),max_score(0);
	int32_t next_anchor(-1);
	for(uint i(start+1);i<template_read.size();++i){
		kmer next(template_read[i]);
		int score(anchors_ordered_according2reads(template_read[start],next,kmer_index));
		if(score>=edge_solidity){
			auto p=longest_ordered_chain_from_anchors(kmer_index,best_chain_computed,i,template_read,edge_solidity);
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

	//~ cerr<<"SCORE of "<<start<<": "<<max_chain+1<<" "<<max_score<<" "<<next_anchor<<endl;
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




vector<kmer> longest_ordered_chain( kmer2localisation& kmer_index,const vector<kmer>& template_read, double edge_solidity){
	unordered_map<uint,score_chain> best_chain_computed;
    vector<kmer> result;
    int32_t max_chain(0),max_score(0);
	int32_t next_anchor(-1);
    for(int32_t i(template_read.size()-1);i>=0;--i){
        auto p=longest_ordered_chain_from_anchors(kmer_index,best_chain_computed,i,template_read,edge_solidity);
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



bool comparable(double x, pair<double,double> deciles){
	//~ return true;
	if(x<deciles.second+5){
		//LOW VALUE
		if(x>deciles.first-5){
			return true;
		}
		if(x/deciles.first<0.5){
			return false;
		}
		return true;
	}else{
		//High value
		if(x/deciles.second>2){
			return false;
		}
		return true;
	}
}



bool comparable(double x,double mean){
	//~ return true;
	if(abs(x-mean)<5){
		return true;
	}
	if(x/mean<0.5 or x/mean>2){
		return false;
	}
	return true;
}



//~ double mean(const vector<uint32_t>& V){
	//~ double first_mean(0);
	//~ uint32_t valid(0);
	//~ for(uint32_t i(0);i<V.size();++i){
		//~ first_mean+=V[i];
	//~ }
	//~ first_mean/=V.size();
	//~ return first_mean;
	//~ double second_mean(0);
	//~ for(uint32_t i(0);i<V.size();++i){
		//~ if(comparable((double)V[i],first_mean)){
			//~ second_mean+=V[i];
			//~ valid++;
		//~ }
	//~ }
	//~ if(valid!=0){
		//~ second_mean/=valid;
	//~ }else{
		//~ return first_mean;
	//~ }
	//~ return second_mean;
//~ }



pair<double,double> deciles( vector<uint32_t>& V){
	sort(V.begin(),V.end());
	return {V[floor((V.size()-1)*0.2)],V[ceil((V.size()-1)*0.8)]};
}



vector<double> average_distance_next_anchor(kmer2localisation& kmer_index,  vector<kmer>& anchors,unordered_map<kmer,uint32_t>& k_count, bool clean){
	vector<double> result;
	vector<uint32_t> v_dis;
	vector<kmer> curated_anchors;
	uint32_t min_distance(5);

	//~ {
		//~ v_dis.clear();
		//~ uint32_t sum(0),count(0);
		//~ auto v_loc1(kmer_index[anchors[i]]);
		//~ uint32_t i1(0);
		//~ while(i1<v_loc1.size()){
			//~ if(v_loc1[i1].read_id==v_loc2[i2].read_id){
				//~ v_dis.push_back(v_loc2[i1].position);
				//~ ++i1;
			//~ }
		//~ }
		//~ auto dec(deciles(v_dis));
		//~ for(uint32_t iD(0);iD<v_dis.size();++iD){
			//~ if(comparable(v_dis[iD],dec)){
				//~ sum+=v_dis[iD];
				//~ count++;
			//~ }
		//~ }
		//~ if(count==0){
			//~ cerr<<"SHOULD NOT HAPPEN"<<endl;
			//~ for(uint32_t iD(0);iD<v_dis.size();++iD){
				//~ sum+=v_dis[iD];
				//~ count++;
			//~ }
		//~ }else{
			//~ result.push_back(sum/count);
		//~ }
	//~ }

	for(uint i(0);i+1<anchors.size();++i){
		v_dis.clear();
		uint32_t sum(0),count(0);
		auto v_loc1(kmer_index[anchors[i]]);//THEY SHOULD BE READS SORTED
		auto v_loc2(kmer_index[anchors[i+1]]);
		//~ sort (v_loc1.begin(), v_loc1.end(), order_according2read_id);
		//~ sort (v_loc2.begin(), v_loc2.end(), order_according2read_id);
		uint32_t i1(0),i2(0);
		while(i1<v_loc1.size() and i2<v_loc2.size()){
			if(v_loc1[i1].read_id==v_loc2[i2].read_id){
				v_dis.push_back(v_loc2[i2].position-v_loc1[i1].position);
				++i1;
				++i2;
			}else if(v_loc1[i1].read_id<v_loc2[i2].read_id){
				i1++;
			}else{
				i2++;
			}
		}
		//~ if(v_dis.empty()){
		//~ }
		auto dec(deciles(v_dis));
		for(uint32_t iD(0);iD<v_dis.size();++iD){
			if(comparable(v_dis[iD],dec)){
				sum+=v_dis[iD];
				count++;
			}
		}
		if(count==0){
			cerr<<"SHOULD NOT HAPPEN"<<endl;
			for(uint32_t iD(0);iD<v_dis.size();++iD){
				sum+=v_dis[iD];
				count++;
			}
		}else{
			result.push_back(sum/count);
		}

		//~ if(count!=0){
			//~ v_sum.push_back(sum);
			//~ v_count.push_back(count);

			//~ result.push_back(sum/count);
		//~ }else{
			//~ cerr<<"SHOULD NOT HAPPEND"<<endl;cin.get();
			//~ result.push_back(-1);
		//~ }
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



vector<vector<string>> split_reads_old(const vector<kmer>& anchors, const vector<double>& relative_positions, const vector<string>& Reads,  kmer2localisation& kmer_index,uint32_t kmer_size){
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



vector<vector<string>> split_reads(const vector<kmer>& anchors, const vector<double>& relative_positions, const vector<string>& Reads,  kmer2localisation& kmer_index,uint32_t kmer_size){
	vector<vector<string>> result(anchors.size()+1);
	if(anchors.size()==0){
		result.push_back(Reads);
		return result;
	}
	for(uint32_t iR(0);iR<Reads.size();++iR){
		//~ cerr<<endl;
		string read=Reads[iR];
		vector<string> split(anchors.size()+1);
		//FIRST AND LAST REGION
		int32_t anchor_position(get_position(kmer_index,anchors[0],iR));
		if(anchor_position!=-1){
			string chunk(read.substr(0,anchor_position));
			//~ if(abs((int)chunk.size()-relative_positions[iA])<get_position(kmer_index,anchors[0],0)*0.5){
			if(comparable(chunk.size(),get_position(kmer_index,anchors[0],0))){
				result[0].push_back(chunk);
			}
		}else{
			//~ result[0].push_back("");
		}
		anchor_position=(get_position(kmer_index,anchors[anchors.size()-1],iR));
		if(anchor_position!=-1){
			string chunk(read.substr(anchor_position));
			if(comparable(chunk.size(),Reads[0].size()-get_position(kmer_index,anchors[anchors.size()-1],0))){
				result[anchors.size()].push_back(chunk);
			}
		}else{
			//~ result[anchors.size()].push_back("");
		}
		for(uint32_t iA(0);iA+1<anchors.size();++iA){
			int32_t anchor_position1(get_position(kmer_index,anchors[iA],iR));
			int32_t anchor_position2(get_position(kmer_index,anchors[iA+1],iR));
			if(anchor_position1!=-1){
				if(anchor_position2!=-1){
					//REGION WITH BOtH ANCHORS
					string chunk(read.substr(anchor_position1,anchor_position2-anchor_position1));
					//~ if(abs((int)chunk.size()-relative_positions[iA])<relative_positions[iA]*0.5){
					if(comparable(chunk.size(), relative_positions[iA])){
						result[iA+1].push_back(chunk);
						//~ cerr<<chunk<<".";
					}else{
						//~ cerr<<"ALIEN"<<endl;
						//~ cerr<<chunk.size()<<" "<<relative_positions[iA]<<endl;
					}
				}else{
					//~ cerr<<'-';
					continue;
					//GOT THE LEFT ANCHOR
					string chunk(read.substr(anchor_position1,relative_positions[iA]));
					if(comparable(chunk.size(),get_position(kmer_index,anchors[0],0))){
						result[iA+1].push_back(chunk);
					}else{
						//~ cerr<<"ALIEN32"<<endl;
					}
				}
			}else{
				if(anchor_position2!=-1){
					//~ cerr<<'-';
					continue;
					//GOT THE RIGHT ANCHOR
					if(anchor_position2>relative_positions[iA]){
						string chunk(read.substr(anchor_position2-relative_positions[iA],relative_positions[iA]));
						if(comparable(chunk.size(),get_position(kmer_index,anchors[0],0))){
							result[iA+1].push_back(chunk);

						}else{
							//~ cerr<<"ALIEN23"<<endl;
						}
					}
				}else{
					//~ cerr<<'-';
					//~ result[iA+1].push_back("");
				}
			}
		}
	}
	return result;
}



int read_string(vector<string>& Vstr,Sequence_T **seq,int do_switch_case,char **comment)
{
	//~ cerr<<"------------READ---------------"<<endl;
  int c,nseq=0,length=0;
  char seq_name[FASTA_NAME_MAX]="",
  line[SEQ_LENGTH_MAX],seq_title[FASTA_NAME_MAX]="";
  char *p;
  stringptr tmp_seq=STRINGPTR_EMPTY_INIT;

  for(uint32_t i(0);i<Vstr.size();++i){

	  if(Vstr[i].empty()){continue;}
	   //~ cerr<<Vstr[i]<<endl;
	char *cstr = new char[Vstr[i].length() + 1];
	strcpy(cstr, Vstr[i].c_str());
	length=Vstr[i].size();
	tmp_seq.p=cstr;

	//TODO DELETE
	if (create_seq(nseq,seq,seq_name,seq_title,tmp_seq.p,do_switch_case))
	  nseq++;
  }
  stringptr_free(&tmp_seq);
  //~ cerr<<nseq<<endl;
  return nseq; /* TOTAL NUMBER OF SEQUENCES CREATED */
}



vector<string> write_string(LPOSequence_T *seq,int nsymbol,char symbol[],int ibundle){
  int i(0);
  vector<string> result;
  int j,nring=0,iprint;
  char **seq_pos=NULL,*p=NULL,*include_in_save=NULL;

  nring=xlate_lpo_to_al(seq,nsymbol,symbol,ibundle, '-',&seq_pos,&p,&include_in_save);
  //~ cerr<<"LPO2al"<<endl;
  LOOPF (i,seq->nsource_seq) { /* NOW WRITE OUT FASTA FORMAT */
    //~ if (ibundle<0 || seq->source_seq[i].bundle_id == ibundle) { /* OR JUST THIS BUNDLE*/
      //~ fprintf(ifile,">%s",seq->source_seq[i].name);
      //~ iprint=0;
      //~ LOOPF (j,nring) { /* WRITE OUT 60 CHARACTER SEQUENCE LINES */
	//~ if (NULL==include_in_save || include_in_save[j]) {
	  //~ fprintf(ifile,"%s%c",iprint%60? "":"\n", seq_pos[i][j]);
	  //~ iprint++; /* KEEP COUNT OF PRINTED CHARACTERS */
	//~ }
		string nadine(seq_pos[i],nring);
		result.push_back(nadine);
		//~ for(uint32_t iN(0);iN<nadine.size();++iN){
			//~ if(nadine[iN]!='-'){
				//~ result.push_back(nadine[iN]);
			//~ }
		//~ }
		//~ break;
      //~ }
      //~ fputc('\n',ifile);
    }

  FREE(p); /* DUMP TEMPORARY MEMORY */
  FREE(include_in_save);
  FREE(seq_pos);
  return result;
}






vector<string> consensus_POA( vector<string>& W){
	 int i,j,ibundle=0,nframe_seq=0,use_reverse_complement=0;
	  int nseq=0,do_switch_case=dont_switch_case,do_analyze_bundles=0;
	  int is_silent = 0;
	  int nseq_in_list=0,n_input_seqs=0,max_input_seqs=100;
	  char score_file[256],seq_file[256],po_list_entry_filename[256],*comment=NULL,*al_name="test align";
	  //~ ResidueScoreMatrix_T score_matrix; /* DEFAULT GAP PENALTIES*/
	  LPOSequence_T *seq=NULL,*lpo_out=NULL,*frame_seq=NULL,*dna_lpo=NULL,*lpo_in=NULL;
	  LPOSequence_T **input_seqs=NULL;
	  FILE *errfile=stderr,*logfile=NULL,*lpo_file_out=NULL,*po_list_file=NULL,*seq_ifile=NULL;
	  char *print_matrix_letters=NULL,*fasta_out=NULL,*po_out=NULL,*matrix_filename=NULL,
		*seq_filename=NULL,*frame_dna_filename=NULL,*po_filename=NULL,*po2_filename=NULL,
		*po_list_filename=NULL, *hbmin=NULL,*numeric_data=NULL,*numeric_data_name="Nmiscall",
		*dna_to_aa=NULL,*pair_score_file=NULL,*aafreq_file=NULL,*termval_file=NULL,
		*bold_seq_name=NULL,*subset_file=NULL,*subset2_file=NULL,*rm_subset_file=NULL,
		*rm_subset2_file=NULL;
	  float bundling_threshold=0.90;
	  int exit_code=0,count_sequence_errors=0,please_print_snps=0,
		report_consensus_seqs=0,report_major_allele=0,use_aggressive_fusion=0;
	  int show_allele_evidence=0,please_collapse_lines=0,keep_all_links=0;
	  int remove_listed_seqs=0,remove_listed_seqs2=0,please_report_similarity;
	  int do_global=0, do_progressive=0, do_preserve_sequence_order=0;
	  char *reference_seq_name="CONSENS%d",*clustal_out=NULL;

  //~ black_flag_init(argv[0],PROGRAM_VERSION);

	matrix_filename="blosum80.mat";
	if(score_matrix_init==false){
		if (read_score_matrix(matrix_filename,&score_matrix)<=0){/* READ MATRIX */
		WARN_MSG(USERR,(ERRTXT,"Error reading matrix file %s.\nExiting", matrix_filename ? matrix_filename: "because none specified"),"$Revision: 1.2.2.9 $");
		}
		score_matrix_init=true;
	}

	//~ cerr<<"GO INSERTION §"<<endl;

	nseq = read_string (W, &seq, do_switch_case, &comment);
	//~ cerr<<"GO INIT AS LPO"<<endl;
	CALLOC (input_seqs, max_input_seqs, LPOSequence_T *);
	for (i=0; i<nseq; i++) {
		//~ cerr<<"i"<<i<<endl;
		input_seqs[n_input_seqs++] = &(seq[i]);
		//~ cerr<<"inputseqnadine"<<endl;
		initialize_seqs_as_lpo(1,&(seq[i]),&score_matrix);//IMPORTANT
		//~ cerr<<"init sucdeees"<<endl;
		if (n_input_seqs == max_input_seqs) {
			max_input_seqs *= 2;
			REALLOC (input_seqs, max_input_seqs, LPOSequence_T *);
		}
	}
	//~ cerr<<"GO CONSENSUS"<<endl;
	lpo_out = buildup_progressive_lpo (n_input_seqs, input_seqs, &score_matrix,use_aggressive_fusion, do_progressive, pair_score_file,matrix_scoring_function, do_global, do_preserve_sequence_order);
	//~ generate_lpo_bundles(lpo_out,bundling_threshold);
	//~ cerr<<"GO OUTPUT"<<endl;
	vector<string> result(write_string(lpo_out,score_matrix.nsymbol,score_matrix.symbol,ibundle));
	//~ cerr<<result<<endl;
	//~ cerr<<"CONSENSUS"<<endl;
	return result;
}


void absoluteMAJ_consensus(vector<string>& V){
	sort(V.begin(),V.end());
	uint score(1),best_occ(0),best_score(0);
	for(uint i(0);i<V.size();++i){
		if(i+1<V.size()){
			if(V[i]!=V[i+1]){
				if(score> 0.5*V.size()){
					V={V[i]};
					return;
				}
				if(score>best_score){
					best_score==score;
					best_occ=i;
				}
				score=1;
			}else{
				score++;

			}
		}else{
			if(score> 0.5*V.size()){
				V={V[i]};
				return;
			}
		}
	}
	//~ V={V[best_occ]};
}



vector<string> easy_consensus(vector<string> V){
	uint32_t non_empty(0);
	//~ absoluteMAJ_consensus(V);
	if(V.size()==1){
		return V;
	}
	//~ uint32_t maximum(0);
	std::set<std::string> mySet;
	for(uint32_t iV(0);iV<V.size();++iV){
		//~ cerr<<V[iV]<<endl;
		//~ maximum=max(maximum,(uint32_t)V[iV].size());
		//~ if(V[iV].size()!=0){
			//~ non_empty++;
			//~ continue;
		//~ }
		mySet.insert(V[iV]);
	}
	// if(V[iV].size()!=V[0].size()){
	if(mySet.size() > 1) {
		V=consensus_POA(V);
		// break;
	} else {
		return {V[0]};
	}
	//~ cerr<<non_empty<<"ne";
	//~ cin.get();
	string result;
	//~ for(uint i(0);i<V.size();++i){
		//~ cerr<<V[i]<<endl;
	//~ }
	//~ cerr<<"END PP"<<endl;
	for(uint32_t iS(0);iS<V[0].size();++iS){
		uint32_t cA,cC,cG,cT,cM;
		cM=cA=cC=cG=cT=0;
		for(uint32_t iV(0);iV<V.size();++iV){
			if(V[iV].size()==0){
				continue;
			}
			switch(V[iV][iS]){
				case 'A': ++cA;break;
				case 'C': ++cC;break;
				case 'G': ++cG;break;
				case 'T': ++cT;break;
				default:
				cM++;
				//~ cerr<<"NOPE"<<V[iV][iS]<<"?"<<endl;
				//~ cerr<<iS<<" "<<V[iV].size()<<" "<<iV<<" "<<V.size()<<endl;
			}
		}
		if(cM>cA and cM>cC and cM>cT and cM>cG){
			// result+=('-');
			continue;
		}
		if(cA>cC and cA>cG and cA>cT){
			result+=('A');
			continue;
		}
		if(cC>cA and cC>cG and cC>cT){
			result+=('C');
			continue;
		}
		if(cG>cA and cG>cC and cG>cT){
			result+=('G');
			continue;
		}
		if(cT>cA and cT>cG and cT>cC){
			result+=('T');
			continue;
		}
		result+=(V[0][iS]);
		// result+='N';
		continue;
		//~ cerr<<"TIE"<<endl;
		return V;
	}
	//~ cerr<<"EASYCONSENSU end"<<endl;

	return {result};
}



vector<vector<string>> global_consensus(const  vector<vector<string>>& V, uint32_t n  ){
	vector<vector<string>> result;
	string stacked_consensus;
	for(uint32_t iV(0);iV<V.size();++iV){
		if(V[iV].size()==0){
			cerr<<"MISSING WINDOWS"<<endl;
			continue;
		}
		vector<string> consensus(easy_consensus(V[iV]));
		//~ cerr<<"EASYCONSENSUS"<<endl;
		//~ cerr<<consensus[0]<<endl;
		//~ if(consensus.size()==1){
			stacked_consensus+=consensus[0];
		//~ }else{
			//~ if(stacked_consensus.size()!=0){
				//~ vector<string> vect(n, stacked_consensus);
				//~ result.push_back(vect);
				//~ stacked_consensus="";
			//~ }
			//~ result.push_back(consensus);
		//~ }
	}
	if(stacked_consensus.size()!=0){
		// stacked_consensus.erase (remove(stacked_consensus.begin(), stacked_consensus.end(), '-'), stacked_consensus.end());
		vector<string> vect(1, stacked_consensus);
		result.push_back(vect);
		stacked_consensus="";
	}
	return result;
}




std::pair<std::vector<std::vector<std::string>>, std::unordered_map<kmer, unsigned>> MSABMAAC(const vector<string>& Reads,uint32_t k, double edge_solidity, unsigned solidThresh){
	int kmer_size(k);
	//~ vector<string> VTest;;
	//~ VTest.push_back("CTGACTGACCCCGTACGTCA");
	//~ VTest.push_back("CTGACTGATTTCGTACGTCA");
	//~ VTest.push_back("CTGACTGAAAACGTACGTCA");
	//~ VTest.push_back("CTGACTGAAAACGTACGTCA");
	//~ VTest.push_back("CTGACTGAAAACGTACGTCAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
	//~ VTest.push_back("CTGACTGATTTCGTACGTCA");
	//~ VTest.push_back("CTGACTGATTTCGTACGTCA");
	//~ VTest.push_back("CTGACTGATTTCGTACGTCA");
	//~ VTest.push_back("CTGACTGATTTCGTACGTCA");
	//~ VTest.push_back("CTGACTGATTTCGTACGTCA");
	//~ VTest.push_back("CTGACTGATTTCGTACGTCA");
	//~ VTest.push_back("CTGACTGCCCCCGTACGTCA");
	//~ VTest.push_back("CTGACTGATTTCGTACGTCA");
	//~ VTest.push_back("CTGACTGATTTCGTACGTCA");
	//~ VTest.push_back("CTGACTGACCCCGTACGTCA");
	//~ auto nadine({VTest});
	//~ auto GC(global_consensus(nadine,1));
	//~ cerr<<GC[0][0]<<endl;
	//~ exit(0);



	kmer2localisation kmer_index;
	std::unordered_map<kmer, unsigned> merCounts;
	fill_index_kmers(Reads,kmer_index,kmer_size,merCounts, solidThresh);
	//~ cerr<<"PHASE 1 done"<<endl;
	//~ return {};

	auto kmer_count(filter_index_kmers(kmer_index,edge_solidity));
	//~ auto kmer_count(filter_index_kmers(kmer_index,percent_shared));
	//~ cerr<<"PHASE 2.1 done"<<endl;
	auto template_read(get_template(kmer_index,Reads[0],kmer_size));
	//~ cerr<<"PHASE 2 done"<<endl;


	vector<kmer> anchors(longest_ordered_chain(kmer_index, template_read,edge_solidity));
	//~ cerr<<"PHASE 3 done"<<endl;


	vector<double> relative_positions=(average_distance_next_anchor(kmer_index,anchors,kmer_count,false));
	//~ cerr<<"PHASE 4 done"<<endl;


	vector<vector<string>> result(split_reads(anchors,relative_positions,Reads,kmer_index,kmer_size));
	//~ cerr<<"PHASE 5 done"<<endl;


	cerr<<""<<result.size()<<"	";
	//~ for(uint i(0);i<result.size();++i){
		//~ for(uint j(0);j<result[i].size();++j){
			//~ cerr<<result[i][j]<<" ";
		//~ }
		//~ cerr<<endl;
	//~ }
	//~ cin.get();

	result=global_consensus(result,Reads.size());
	//~ cerr<<"PHASE 6 done"<<endl;

	return std::make_pair(result, merCounts);
}
