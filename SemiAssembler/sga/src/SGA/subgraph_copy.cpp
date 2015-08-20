//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// subgraph - extract a subgraph from an assembly graph
//
#include <iostream>
#include <cstring>
#include <sstream>
#include <vector>
#include <queue>
#include <fstream>
#include <stdio.h> 
#include <omp.h>
#include "subgraph.h"
#include "Util.h"
#include "SGUtil.h"
#include "SGAlgorithms.h"
#include "SGVisitors.h"
#include "Timer.h"
#include <time.h>
//================================================== 
#include "SuffixArray.h"
#include "SGACommon.h"
#include "BWTAlgorithms.h"
#include "ReadInfoTable.h"
#include "BWT.h"
#include "BWTIntervalCache.h"
#include "OverlapCommon.h"
#include "SGSearch.h"
#include "Bigraph.h"

//================================================== 
//typedef std::map<VertexID, Vertex*> VertexPtrMap;
//typedef __gnu_cxx::hash_map<VertexID, Vertex*> VertexPtrMap;
//typedef std::tr1::unordered_map<VertexID, Vertex*> VertexPtrMap;
#include <tr1/unordered_map>
#include <tr1/unordered_set>
#define HashMap std::tr1::unordered_map
#define HashSet std::tr1::unordered_set
typedef std::tr1::hash<std::string> StringHasher;
#include <sparsehash/dense_hash_set>
#define DenseHashSet google::dense_hash_set
//typedef DenseHashSet<std::string, StringHasher> EdgeHashSet;
typedef DenseHashSet<VertexID, StringHasher> EdgeHashSet;
typedef EdgeHashSet::iterator EdgeHashSetIter;
//==================================================
#define OMP 63

// functions
void addNeighborsToSubgraph_dfs(Vertex* pCurrVertex, StringGraph* pSubgraph, EdgeDir dir, int *stop, int *Turn_off_flag, unsigned int *final_distance, unsigned int upper_distance, unsigned int curr_distance, unsigned int *leaf_count, unsigned int *VertexCount, std::string &LargeDegreeVertex, EdgeHashSet* End_set, EdgeHashSet* visit_Vertex_set);
void load_MpSeqInfo();  // load matepair file
void OutputContig(StringGraph* pSubgraph, int Mode_Control, unsigned int mateNum, unsigned int final_distance, int thread_num, std::string outStartVertexID, std::string outEndVertixID, int Comp_flagS, int Comp_flagE);
//void OutputContig_bfs(StringGraph* pSubgraph,int Case_Control,unsigned int mateNum,unsigned int total_distance,int thread_num,std::string outStartVertexID,std::string outEndVertexID);
void OutputContig_bfs_case1(StringGraph* pSubgraph,unsigned int mateNum,unsigned int total_distance,int thread_num,std::string outStartVertexID);//for case 1
void OutputContig_bfs_case2(StringGraph* pSubgraph_start,StringGraph* pSubgraph_end,unsigned int mateNum,int status,int thread_num,std::string outStartVertexID,std::string outEndVertexID);//for case 2
void copyVertexToSubgraph(StringGraph* pSubgraph, const Vertex* pVertex);
void addNeighborsToSubgraph(Vertex* pCurrVertex,StringGraph* pSubgraph,int span);
std::string int2str(int &); // int to string

// Getopt
#define SUBPROGRAM "subgraph"
static const char *SUBGRAPH_VERSION_MESSAGE =
SUBPROGRAM " Version " PACKAGE_VERSION "\n"
"Written by Yu-Han Su.\n"
"\n"
"\n";

static const char *SUBGRAPH_USAGE_MESSAGE =
"Usage: " PACKAGE_NAME " " SUBPROGRAM " [OPTION] ... PE_index_prefix Mp_seq.fq ASQGFILE\n"
"Extract the subgraph around the sequence with ID from an asqg file with direction.\n"
"\n"
"  -v, --verbose                        display verbose output\n"
"      --help                           display this help and exit\n"
"      -o, --out=FILE                   write the subgraph to FILE (default: subgraph.asqg.gz, subgraph2.asqg.gz)\n"
"      -s, --size=N                     the size of the subgraph to extract, all vertices that are at most N hops\n"
"                                       away from the root will be included (default: 5)\n"
"\n"
"Search Control parameters:\n"
"      -d                               0: sense, 1: antisense, 2: bidirectional\n"
"      -m, --size=N                     mean of insert size (default: 3000)\n"
"      -a, --size=N                     variance of insert size (default: 1000)\n"
"      -t, --size=N                     search depth (default: mean + variance)\n"
"\n"
"      -F, --size=N                     Upper bound of Leaf count   (default: 0(all) )\n"
"      -V, --size=N                     Upper bound of Vertex count (default: 0(all) )\n"
"\n"
"Trimming parameters:\n"
"      -X, --cut-terminal=N             cut off terminal branches in N rounds (default: 10)\n"
"      -L, --min-branch-length=LEN      remove terminal branches only if they are less than LEN bases in length (default: 150)\n"
"          --transitive-reduction       remove transitive edges from the graph. Off by default.\n"
"          --fastq						input matepair (fastq format).\n"
"          --fasta						input matepair (fasta format).\n"
"\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

namespace opt
{
	static unsigned int verbose;
	static std::string asqgFile;            // Input asqg
	static std::string indexPrefix;         // Input pair-end index prefix
	static std::string MpFile;	            // Input mate-pair file

	static std::string prefix;               // Output prefix
	static std::string outSubgraph;          // Output Subgraph
	static std::string outContigsFile;       // Output contigs file

	// search path parameters
	static unsigned int Edir = 0;              // 0: sense, 1: antisense, 2: bidirectional
	static unsigned int span = 5;              // subgrph depth
//	static unsigned int mean_d = 4500;         // mean of insert size
//	static unsigned int var_d = 1500;          // variance of insert size
	 static unsigned int mean_d = 1900;         // mean of insert size
        static unsigned int var_d = 100;
	// Search Control parameters
	//static int depth = mean_d + var_d + 200;   // subgraph depth
	static int depth = 1000000;
	static unsigned int LeafCount = 0;
	static unsigned int TotalVertexCount = 0;

	// Trim parameters
	static int numTrimRounds = 5;
	static unsigned int trimLengthThreshold = 300;

	// Vistor Control
	static bool bPerformTR  = false;            // SGTransitiveReductionVisitor
	static bool loadFastQ   = false;
	static bool loadFastA   = false;

	// Data array
	static std::vector<std::string> MpInfo1_ID;     // mp_ID1_array
	static std::vector<std::string> MpInfo1_SEQ;    // mp_ID1_array
	static std::vector<std::string> MpInfo2_ID;     // mp_ID2_array
	static std::vector<std::string> MpInfo2_SEQ;    // mp_ID2_array
        static std::vector<std::string> MpInfo_length;
}

static const char* shortopts = "o:s:d:m:a:t:F:L:X:V:";

enum { OPT_HELP = 1, OPT_VERSION, OPT_TR, OPT_FA, OPT_FQ };

static const struct option longopts[] = {
	{ "verbose",                  no_argument,       NULL, 'v' },
	{ "prefix",                   required_argument, NULL, 'o' },
	{ "span",                     required_argument, NULL, 's' },  // span
	{ "Edir",                     required_argument, NULL, 'd' },  // Edir
	{ "mean-of-insert-size",      required_argument, NULL, 'm' },  // mean of insert size
	{ "variance-of-insert-size",  required_argument, NULL, 'a' },  // variance of insert size
	{ "search-depth",             required_argument, NULL, 't' },  // search depth
	{ "leaf-count",               required_argument, NULL, 'F' },  // leaf_count
	{ "vertex-count",             required_argument, NULL, 'V' },  // leaf_count
	{ "cut-terminal",             required_argument, NULL, 'X' },  // numTrimRounds
	{ "min-branch-length",        required_argument, NULL, 'L' },  // trimLengthThreshold
	{ "help",                     no_argument,       NULL, OPT_HELP },
	{ "version",                  no_argument,       NULL, OPT_VERSION },
	{ "transitive-reduction",     no_argument,       NULL, OPT_TR },
	{ "fastq",                    no_argument,       NULL, OPT_FQ },
	{ "fasta",                    no_argument,       NULL, OPT_FA },
	{ NULL, 0, NULL, 0 }
};

//
// Main
//
int subgraphMain(int argc, char** argv)
{
	Timer* pTimer = new Timer("Total Time");
	parseSubgraphOptions(argc, argv);
	load_MpSeqInfo();
	subgraph();
	delete pTimer;

	return 0; 
}

void subgraph()
{
        std::fstream insertion_file;
        insertion_file.open ("Assemble_inserions.fa",std::ios::out);
      //  std::ifstream insertion_file("Assemble_inserions.fa",std::ios::out);
        if(!insertion_file)

        {

                std::cerr << "Can't open file!\n";

                exit(1);

        }
        
	//
	// Visitor functors
	SGGraphStatsVisitor statsVisit;
	SGContainRemoveVisitor containVisit;
	SGTransitiveReductionVisitor trVisit;
	SGTrimVisitor trimVisit(opt::trimLengthThreshold);

	//
	// load ASQG file (matepair map paired reads)
	clock_t t_start, t_end;
	t_start = clock();
	StringGraph* pGraph = SGUtil::loadASQG(opt::asqgFile, 0, true);
	t_end = clock();
	std::cout << "Done: loadASQG " << "\tTime:" << ((t_end-t_start)/(double)(CLOCKS_PER_SEC)) << std::endl;
	pGraph->visit(statsVisit);
	//pGraph->writeASQG("newasqg_1.asqg");
	//pGraph->writeASQG();
	//
	// Remove containments from the graph
	t_start = clock();
	while(pGraph->hasContainment())
		pGraph->visit(containVisit);
	t_end = clock();
	std::cout << "Done: SGContainRemoveVisitor" << "\tTime:" << ((t_end-t_start)/(double)(CLOCKS_PER_SEC)) << std::endl;
	pGraph->visit(statsVisit); 

	//
	// Remove any extraneous transitive edges that may remain in the graph
	if(opt::bPerformTR)
	{
		t_start = clock();
		pGraph->visit(trVisit);
		t_end = clock();
		std::cout << "Done: SGTransitiveReductionVisitor" << "\tTime:" << ((t_end-t_start)/(double)(CLOCKS_PER_SEC)) << std::endl;
		pGraph->visit(statsVisit); 
	}
//pGraph->writeASQG("newasqg.asqg");
	//
	// load FM index and reads file (for search MP reads mapping PE reads)
	t_start = clock();
	int sampleRate = 128;	//sga default
	BWT* pBWT = new BWT(opt::indexPrefix + BWT_EXT, sampleRate);	// load *.bwt
	BWT* pRBWT = NULL;
	SampledSuffixArray* pSSA = NULL;

	int intervalCacheLength=10;	//sga default
	BWTIntervalCache* pIntervalCache = new BWTIntervalCache(intervalCacheLength, pBWT);

	BWTIndexSet indexSet;
	indexSet.pBWT = pBWT;
	indexSet.pRBWT = pRBWT;
	indexSet.pSSA = pSSA;
	indexSet.pCache = pIntervalCache;

	SampledSuffixArray* pSSA2 = new SampledSuffixArray(opt::indexPrefix + SAI_EXT, SSA_FT_SAI);	//load *.sai
	std::string input_pe_file = opt::indexPrefix + ".fa" ; // sga_pe_rmdup.fa
	ReadInfoTable* pTargetReadInfoTable = new ReadInfoTable(input_pe_file);

	t_end = clock();
	std::cout << "Done: Load FM index" << "\tTime:" << ((t_end-t_start)/(double)(CLOCKS_PER_SEC)) << std::endl;

	// matepair path visitor (for each matepair)
	unsigned int mateNum=0;
	unsigned int mateNum_up_bound=opt::MpInfo1_ID.size();
	//===========================================================
	//#pragma omp parallel for num_threads(30) private(mateNum) //開thread
	//===========================================================
	for(mateNum=0 ; mateNum < mateNum_up_bound ; mateNum++)  // for each mate pair
	{	
		clock_t m_start,m_end;	//Timer* pTimerE = new Timer("Each matepair"); 
		m_start = clock();

		std::string StartVertixID = opt::MpInfo1_ID[mateNum];  // start ID in mp_ID1 (xxx/1)
		std::string EndVertixID   = opt::MpInfo2_ID[mateNum];  // end   ID in mp_ID2 (xxx/2)

		std::vector<std::string> StartVertixID_list;
		std::vector<std::string> EndVertixID_list;
		std::vector<std::string> LDFS;
		std::vector<std::string> LDFE;
		
		int flag_end; //0 means sense ,1 means antisense
		int flag_start;

		//===========================================================================================
		//
		// Get EndVertexID and insert to EndVetex_hash
		// get interval(forward)
		int substr_len = opt::MpInfo2_SEQ[mateNum].length();
		BWTInterval interval,interval_rev;
		std::string kmer;
		
                std::string predict_lenS= opt::MpInfo_length[mateNum];
                int predict_len = std::atoi(predict_lenS.c_str());
                std::cout << "predict_len= " << predict_len << "\n";		
		//for(int cut_site=0; cut_site<substr_len-25; cut_site+=3)//old way
		for(int cut_size=substr_len;cut_size>=90;cut_size-=3)//new way
		{
			//for(int cut_len=1; cut_len < substr_len-25; cut_len+=3)
			for(int str_start=0;str_start<substr_len-cut_size+1;str_start+=3)
			{
				//kmer = opt::MpInfo2_SEQ[mateNum].substr(cut_site,substr_len-cut_len);   // EndVertex_seq
				kmer = opt::MpInfo2_SEQ[mateNum].substr(str_start,cut_size);
				interval = BWTAlgorithms::findInterval(indexSet, kmer);	 // Find FM index interval for EndVertex_SEQ
				//std::cout << "mp" << mateNum+1 << " cut:" << cut_site <<" End_kmer = " << kmer << "\tlen=" << kmer.length() << " find:" << interval.size() << std::endl;
				if (interval.size() > 0 )
					break;
			}
			if (interval.size() > 0 )
				break;
		}
	
		//std::cout<<"END interval size= "<<interval.size()<<std::endl;
		// get interval(reverse)
		
		std::string MpInfo2_seq_rev;
		MpInfo2_seq_rev = reverseComplement(opt::MpInfo2_SEQ[mateNum].c_str());	// reverse
		//for(int cut_site=0; cut_site<substr_len-25; cut_site+=3)
		for(int cut_size=substr_len;cut_size>=90;cut_size-=3)
		{
			//for(int cut_len=1; cut_len < substr_len-25; cut_len+=3)
			for(int str_start=0;str_start<substr_len-cut_size+1;str_start+=3)
			{
				//kmer = MpInfo2_seq_rev.substr(cut_site,substr_len-cut_len);		// EndVertex_seq
				//kmer = opt::MpInfo2_seq_rev[mateNum].substr(str_start,cut_size);
				kmer = MpInfo2_seq_rev.substr(str_start,cut_size);
				interval_rev = BWTAlgorithms::findInterval(indexSet, kmer);		// Find FM index interval for EndVertex_SEQ_reverse
				//std::cout <<"cut_len=" <<cut_len<<"mp" << mateNum+1 << " cut:" << cut_site <<" End_kmer = " << kmer << "\tlen=" << kmer.length() << " find:" << interval_rev.size() << std::endl;
				if (interval_rev.size() > 0 )
					break;
			}
				if (interval_rev.size() > 0 )
					break;
		}
	//	std::cout<<"END rev interval size= "<<interval_rev.size()<<std::endl;
		//std::cout<<"mp ="<<mateNum+1<<"is out\n";
		 
		//===========================================================================================
		//
		// Get EndVertexID then insert to EndVetex_hash
		EdgeHashSet* End_set = new EdgeHashSet;  // EdgeHash for EndVertex (Start to End)
		End_set->set_empty_key("");
		flag_end=0;
		if( ( interval.size() + interval_rev.size() ) > 0 )		// MP hit num != 0
		//if( interval.size() > 0)		// MP hit num != 0
		{
			int count=0;
			if(interval.size()==0)//mean interval_rev.size()>0
				flag_end=1; 
			for (int64_t j=interval.lower ; j<=interval.upper ; ++j)
			{
				if(count>5) break; // only keep 5 readID
				int64_t saIdx = j;
				const ReadInfo& targetInfo = pTargetReadInfoTable->getReadInfo(pSSA2->calcSA(saIdx,pBWT).getID()); // get PE id
				EndVertixID=targetInfo.id;
				//std::cout<<"END="<<EndVertixID<<std::endl;
				End_set->insert(EndVertixID.c_str());    // insert to End_point Hash_set
				EndVertixID_list.push_back(EndVertixID); // push to EndVertex_list
				//std::cout << count << "\tEndVertixID = " << EndVertixID << std::endl;
				count++;
			}
			
			 // count=0;
			for (int64_t j=interval_rev.lower ; j<=interval_rev.upper ; ++j)
			   {
			   if(count>5) break; // only keep 5 readID
			   int64_t saIdx = j;
			   const ReadInfo& targetInfo = pTargetReadInfoTable->getReadInfo(pSSA2->calcSA(saIdx,pBWT).getID());
			   EndVertixID=targetInfo.id;
			  // std::cout<<"END="<<EndVertixID<<std::endl;
			   End_set->insert(EndVertixID.c_str());    // End_point Hash_set
			   EndVertixID_list.push_back(EndVertixID);
			//std::cout << count << "\trevEndVertixID = " << EndVertixID << std::endl;
			count++;
			} 
			 
		}
		
		else 
		{
			std::cout << "mp" << mateNum+1 << "\t End Vertex not found.\n";
			continue;
		}
		//===========================================================================================
		//
		// Get StartVertexID list 
		// get interval(forward)
		substr_len = opt::MpInfo1_SEQ[mateNum].length();
		//std::cout << "FUCK"<<std::endl;
		//for(int cut_site=0; cut_site<substr_len-50; cut_site+=3)
		for(int cut_size=substr_len;cut_size>=90;cut_size-=3)
		{
			//for(int cut_len=1; cut_len < substr_len-25; cut_len+=3)
			for(int str_start=0;str_start<substr_len-cut_size+1;str_start+=3)
			{
				//kmer = opt::MpInfo1_SEQ[mateNum].substr(cut_site,substr_len-cut_len);
				//kmer = opt::MpInfo1_SEQ[mateNum].substr(cut_site,substr_len);// for CYT
				kmer = opt::MpInfo1_SEQ[mateNum].substr(str_start,cut_size);
				interval = BWTAlgorithms::findInterval(indexSet, kmer);	 // Find FM index interval for StartVertex_SEQ (kmer)
				//std::cout << "mp" << mateNum+1 << " cut:" << cut_site <<" Start_kmer = " << kmer << "\tlen=" << kmer.length() << " find:" << interval.size() << std::endl;
				if (interval.size() > 0 )
					break;
			}
			if (interval.size() > 0 )
				break;
		}
		//std::cout<<"Start interval size= "<<interval.size()<<std::endl;
		// get interval(reverse)
		
	
		std::string MpInfo1_seq_rev;
		MpInfo1_seq_rev = reverseComplement(opt::MpInfo1_SEQ[mateNum].c_str());	
		//for(int cut_site=0; cut_site<substr_len-50; cut_site+=3)
		for(int cut_size=substr_len;cut_size>=90;cut_size-=3)
		{
			//for(int cut_len=1; cut_len < substr_len-25; cut_len+=3)
			for(int str_start=0;str_start<substr_len-cut_size+1;str_start+=3)
				{
					//kmer = MpInfo2_seq_rev.substr(cut_site,substr_len-cut_len);		// EndVertex_seq
					//kmer = MpInfo2_seq_rev.substr(cut_site,substr_len);//for CYT
					kmer = MpInfo1_seq_rev.substr(str_start,cut_size);
					interval_rev = BWTAlgorithms::findInterval(indexSet, kmer);		// Find FM index interval for EndVertex_SEQ_reverse
					//std::cout << "mp" << mateNum+1 << " cut:" << cut_site <<" End_kmer = " << kmer << "\tlen=" << kmer.length() << " find:" << interval_rev.size() << std::endl;
					if (interval_rev.size() > 0 )
						break;
				}
					if (interval_rev.size() > 0 )
						break;
		}
		// std::cout<<"Start interval rev size= "<<interval_rev.size()<<std::endl;
		//===========================================================================================
		// Get StartVertexID list then insert to StartVertixID_list
		int SelfVertex=0;
		flag_start=0;
		if( ( interval.size() + interval_rev.size() ) > 0 )		// MP hit num != 0
		//if( interval.size() > 0)		// MP hit num != 0
		{
			int count=0;
			if(interval.size()==0)//mean interval_rev.size()>0
				flag_start=1;
			for (int64_t j=interval.lower ; j<=interval.upper ; ++j)
			{
				if(count>5) break;	// only keep 5 readID
				int64_t saIdx = j;
				const ReadInfo& targetInfo = pTargetReadInfoTable->getReadInfo(pSSA2->calcSA(saIdx,pBWT).getID()); // get PE id
				StartVertixID=targetInfo.id;
				StartVertixID_list.push_back(StartVertixID);	// push to StartVertex_list
				//std::cout << count << "\tStartVertixID = " << StartVertixID << std::endl;
//dont know
				//EdgeHashSetIter RHSiter = End_set->find(StartVertixID);		// check StartVertex in EndVertex set.
			//	if (RHSiter != End_set->end())   // StartVertex == EndVertex
				//	SelfVertex=1;
//dont know
				count++;
			}
				//count=0;
				for (int64_t j=interval_rev.lower ; j<=interval_rev.upper ; ++j)
				{
				if(count>5) break;	// only keep 5 readID
				int64_t saIdx = j;
				const ReadInfo& targetInfo = pTargetReadInfoTable->getReadInfo(pSSA2->calcSA(saIdx,pBWT).getID());
				StartVertixID=targetInfo.id;
				StartVertixID_list.push_back(StartVertixID);
				//EdgeHashSetIter RHSiter = End_set->find(StartVertixID);		// check StartVertex in EndVertex set.
				//if (RHSiter != End_set->end())   // StartVertex == EndVertex
					//SelfVertex=1;
					
				count++;
				}
		}
		else 
		{
			std::cout << "mp" << mateNum+1 << "\t Start Vertex not found.\n";
			continue;
		}
		/*Choose the StartVertex and EndVertex*/
		int same_flag=0;
		
		for(size_t i=0 ;i<EndVertixID_list.size();i++)
		{
			EndVertixID=EndVertixID_list[i];
			for(size_t j=0;j<StartVertixID_list.size();j++)
			{
				StartVertixID=StartVertixID_list[j];
				if(EndVertixID!=StartVertixID)
				{
					same_flag=1;
					break;
				}	
			}
			if(same_flag==1)
				break;
		
		}
		
		if(same_flag==0)//if same then give it up
			{
				std::cout << "mp" << mateNum+1 << "\t StartVertex == EndVertex.\n";
                                std::cout << "contig id= "<< StartVertixID <<std::endl;
				std::cout << pGraph->getVertex(StartVertixID)->getStr() << std::endl;
				insertion_file << opt::MpInfo1_ID[mateNum] << "\n" << pGraph->getVertex(StartVertixID)->getStr() << "\n";
				continue;
			}
		StringGraph* pSubgraph_same = new StringGraph;
	/*	Vertex* sameVertex = NULL;
		sameVertex=pGraph->getVertex(StartVertixID);
		if(SelfVertex==1)
		{
			std::cout << "mp" << mateNum+1 << "\t StartVertex == EndVertex.\n";
			//std::cout<<"Start="<<StartVertixID<<"End="<<EndVertixID<<std::endl;
			//std::cout<<"mp "<<mateNum+1<<std::endl;
			//std::cout<<sameVertex->getStr()<<std::endl;
			continue;
		}*/
		//std::cout << "mp = " <<mateNum+1 <<"\tfuck"<< std::endl;

		// Set the graph parameters to match the main graph
		StringGraph* pSubgraph = new StringGraph;
		pSubgraph->setContainmentFlag(pGraph->hasContainment());
		pSubgraph->setTransitiveFlag(pGraph->hasTransitive());
		pSubgraph->setMinOverlap(pGraph->getMinOverlap());
		pSubgraph->setErrorRate(pGraph->getErrorRate());

		// Record visit_Vertex
		EdgeHashSet* visit_Vertex_set = new EdgeHashSet;  // DFS record visted vertex
		visit_Vertex_set->set_empty_key("");

		// initail search control parameter
		Vertex* pRootVertex = NULL;
		int stop=0;
		unsigned int final_distance=0, upper_distance=0, curr_distance=0;
		unsigned int leaf_count=0;
		unsigned int VertexCount=0;
		std::string LargeDegreeVertexFromS="";
		std::string LargeDegreeVertexFromE="";
		std::string outStartVertexID="",outEndVertexID="";
		int Mode_Control=0;   // mode0: SE successful, mode2 = SE failed but ES sucessful, mode3 = SE & ES failed (convert to Nver SuperReads)
		int Comp_flagS=-1,Comp_flagE=-1;	// Reads Direction Check
		int Turn_off_flag=0; // flag to terminate DFS

		unsigned int VertexCount1=0;
		unsigned int VertexCount2=0;
		unsigned int VertexCount3=0;
		unsigned int VertexCount4=0;

		/* Part II. using BFS to find the path from StartVertex to EndVertex */
		
		//std::cout<<"Start="<<StartVertixID<<" End="<<EndVertixID<<std::endl;
		int thread_num=omp_get_thread_num();
		Vertex* pX=pGraph->getVertex(StartVertixID);
		Vertex* pY=pGraph->getVertex(EndVertixID);
		std::cout << pGraph->getVertex(StartVertixID)->getStr() << std::endl;
                std::cout << pGraph->getVertex(EndVertixID)->getStr() << std::endl;

		SGWalkVector allwalk;
		SGWalkVector StartWalks;
		SGWalkVector EndWalks;
	//	StringGraph* pSubgraph_start = new StringGraph;
	//	StringGraph* pSubgraph_end = new StringGraph;
	//	StringGraph* pSubgraph_all = new StringGraph;
		//std::cout<<"sense="<<flag_start<<" pX len="<<pX->getSeqLen()<<" pY len="<<pY->getSeqLen()<<std::endl;
		if(flag_start==0)
			SGSearch::findWalks(pX,pY,ED_SENSE,predict_len+opt::var_d+pX->getSeqLen()+pY->getSeqLen(),opt::LeafCount,false,allwalk);//change 500->opt::LeafCount
		//	SGSearch::findWalks(pX,pY,ED_SENSE,opt::mean_d+opt::var_d,opt::LeafCount,false,allwalk);
		else if(flag_start==1)    
            	        SGSearch::findWalks(pX,pY,ED_ANTISENSE,predict_len+opt::var_d+pX->getSeqLen()+pY->getSeqLen(),opt::LeafCount,false,allwalk);//
		//	SGSearch::findWalks(pX,pY,ED_ANTISENSE,opt::mean_d+opt::var_d,opt::LeafCount,false,allwalk);
		int flag=0;
		int status=0;
		int total_distance;
		std::cout<<"allwalk="<<allwalk.size()<<" pX= "<<StartVertixID<<" pY="<<EndVertixID<<std::endl;
		//std::cout<<"Sense="<<ED_SENSE<<std::endl;
		//std::cout<<"ANTISense="<<ED_ANTISENSE<<std::endl;
		std::string walkStr="";
		int goal_path=5000;
		int goal_path_walk;
                int goal_path_cmp;
		if(allwalk.size()>0)
		{
			/* need to choose the path
			*/
			size_t temp;
			int flag_path;//check there exists a legal path between Start and End or not.
			flag_path=0;
			std::cout << "mean=" << predict_len << "\t var=" << opt::var_d << std::endl ;
			for(size_t i=0;i<allwalk.size();i++)// not sure
			{
				//if((allwalk[i].getString(SGWT_INTERNAL).size())>=(opt::mean_d-opt::var_d)&&(allwalk[i].getString(SGWT_INTERNAL).size())<=(opt::mean_d+opt::var_d))
			/*	if(allwalk[i].getString(SGWT_INTERNAL).size()<=(opt::mean_d+opt::var_d))
				{	temp=i;
					flag_path=1;
					break;
				}*/
				
				//if((allwalk[i].getString(SGWT_START_TO_END).size()>(opt::mean_d+opt::var_d+pX->getSeqLen()+pY->getSeqLen()))||(allwalk[i].getString(SGWT_START_TO_END).size()<(opt::mean_d-opt::var_d+pX->getSeqLen()+pY->getSeqLen())))//make sure the path is legal or not
				//if((allwalk[i].getString(SGWT_INTERNAL)>(opt::mean_d+opt::var_d))||(allwalk[i].getString(SGWT_INTERNAL)<(opt::mean_d-opt::var_d)))
				if(!(allwalk[i].getString(SGWT_INTERNAL).size()<=(predict_len+opt::var_d)&&(allwalk[i].getString(SGWT_INTERNAL).size()>=(predict_len-opt::var_d))))
				{
				//	allwalk.erase(allwalk.begin()+i);// not sure
					std::cout << "allwalk-" << i << "=" << "\t illegal path. Size="<< allwalk[i].getString(SGWT_INTERNAL).size() << std::endl ;
				}
				else
				{
					std::cout << "allwalk-" << i << "=" << "\t legal path.  Size="<< allwalk[i].getString(SGWT_INTERNAL).size() << std::endl ;
				}
                                goal_path_cmp=(allwalk[i].getString(SGWT_INTERNAL).size())-(predict_len);
				if(goal_path_cmp < 0) goal_path_cmp = -goal_path_cmp;
				if(goal_path_cmp < goal_path)
				{
                                        goal_path=goal_path_cmp;
					goal_path_walk=i;
				}
	
		
			}
                        std::cout <<"goal_path_walk= " << goal_path_walk << std::endl ;
			if(allwalk.size()==0)
			{
				//std::cout<<"no path"<<std::endl;
				//std::cout<<"Illigal size="<<allwalk[0].getString(SGWT_INTERNAL).size()<<std::endl;
				std::cout << "mp" << mateNum+1 << "_case1" << "\t illegal path." << std::endl ;
				continue;
			}
	//		VertexPtrVec vertexPtrVector_all;
			//for(size_t i=0;i<allwalk.size();i++)
	//		for(size_t k=0;k<allwalk.size();k++)
	//		{
			for(size_t k=goal_path_walk;k<=goal_path_walk;k++)
			{
				
				//adddddddddd
				std::cout << "w=" << k << std::endl ;
				SGWalk& walk = allwalk[k];
				std::string str = walk.getString(SGWT_START_TO_END);
				std::cout << "Str: " << str << "\n";
				walk.print();
				std::cout << "wwwww" << std::endl ;
				insertion_file << opt::MpInfo1_ID[mateNum] << "\n" << str << "\n";
				//addddddddd

				
				std::cout << "k1=" << k << std::endl ;
                		StringGraph* pSubgraph_all = new StringGraph;
				VertexPtrVec vertexPtrVector_all;
				VertexPtrVec vertex_tmp=allwalk[k].getVertices();
				//vertexPtrVector_all.push_back(vertex_tmp);
				vertexPtrVector_all.insert(vertexPtrVector_all.end(),vertex_tmp.begin(),vertex_tmp.end());
				for(size_t i=0;i<vertexPtrVector_all.size();i++)
                       	 	{	
                                	copyVertexToSubgraph(pSubgraph_all,vertexPtrVector_all[i]);
					
				

                        	}
			
                        	for(size_t i=0;i<vertexPtrVector_all.size()-1;i++)
                        	{
                                	EdgePtrVec edges=vertexPtrVector_all[i]->getEdges();
                                	for(size_t j=0;j<edges.size();j++)
                                	{
                                        	if(edges[j]->getEnd()==vertexPtrVector_all[i+1])
                                        	{
                                	                Overlap ovr = edges[j]->getOverlap();
                        	                        SGAlgorithms::createEdgesFromOverlap(pSubgraph_all,ovr,true);
							//addd

                	                        }
        	                        }
					
	                        }
			
				SGDuplicateVisitor dupVisit;
			//	pSubgraph_all->visit(dupVisit);
				flag=1;
				if(k<20 && k > 10){
                                char char_mpNum[20];
                               // sprintf(char_mpNum,"%d",mateNum+1);
				sprintf(char_mpNum,"%d",k);
                                std::string string_mpNum = char_mpNum ;
                                std::string ContigName = "mp_" + string_mpNum + ".dot";
                               
                                pSubgraph_all->writeDot(ContigName);}
				std::cout << "k2=" << k << std::endl ;
				OutputContig_bfs_case1(pSubgraph_all,mateNum,total_distance,thread_num,StartVertixID);
				delete pSubgraph_all;
				std::cout << "k3=" << k << std::endl ;

			}
		//	std::sort(vertexPtrVector_all.begin(),vertexPtrVector_all.end());
			//std::vector<int>::iterator it;
			//std::unique(vertexPtrVector_all.begin(),vertexPtrVector_all.end());
			//vertexPtrVector_all.resize(std::distance(vertexPtrVector_all.begin(),it));
		/*	for(size_t i=0;i<vertexPtrVector_all.size();i++)
			{
				copyVertexToSubgraph(pSubgraph_all,vertexPtrVector_all[i]);
			}
			for(size_t i=0;i<vertexPtrVector_all.size()-1;i++)
			{
				EdgePtrVec edges=vertexPtrVector_all[i]->getEdges();
				for(size_t j=0;j<edges.size();j++)
				{
					if(edges[j]->getEnd()==vertexPtrVector_all[i+1])
					{
						Overlap ovr = edges[j]->getOverlap();
						SGAlgorithms::createEdgesFromOverlap(pSubgraph_all,ovr,true);
					}
				}
			}
			
	
				flag=1;
				char char_mpNum[20];
				sprintf(char_mpNum,"%d",mateNum+1);
				std::string string_mpNum = char_mpNum ; 
				std::string ContigName = "mp_" + string_mpNum + ".dot";
				pSubgraph_all->writeDot(ContigName);
		*/		
		}
		
		
		else //can't find
		{
			std::cout << "mp" << mateNum+1 << "_case1" << "\t Can not find." << std::endl ;
		}
		/* Part III. Scaffold */
		if(flag==1 || same_flag==0)//case 1 S--->E
		{//(StringGraph* pSubgraph,int Case_Control,unsigned int mateNum,unsigned int total_distance,int thread_num,std::string outStartVertexID,std::string outEndVertexID);
			
		//	OutputContig_bfs_case1(pSubgraph_all,mateNum,total_distance,thread_num,StartVertixID);
		}
		
	/*	delete pSubgraph_all;
		delete pSubgraph_start;
		delete pSubgraph_end;
	*/	//myBFS
		
		
	/*	
		
		for(size_t i = 0; i < StartVertixID_list.size(); i++)    // for each start point // Start to End search
		{
			//std::cout<<"i="<<i<<"\tStart="<<StartVertixID_list[i]<<std::endl;
			//std::cout << "SEmod:::SEmode Step.1 " << std::endl;
			if(i>10) break;	// only try three StartVertex
			StartVertixID=StartVertixID_list[i];	// Get the start(root) vertex from main graph
			//std::cout<<"mp="<<mateNum+1<<"\t Start= "<<StartVertixID<<".size="<<StartVertixID_list.size()-1<<"\n";
			pRootVertex = pGraph->getVertex(StartVertixID);
			curr_distance=0, final_distance=0;
			stop=0, Turn_off_flag=0;
			leaf_count=0, VertexCount=0;
			LargeDegreeVertexFromS = "";	// no necessary initial because only keep one time
			upper_distance=opt::depth; // max of search depth

			if(pRootVertex == NULL)
			{
				std::cout << "mp" << mateNum+1 << "\t Step.1 StartVertex:: " << StartVertixID << " not found in the graph.\n";
				continue;
			}
			else
			{
				//EdgeDir StartVertixID_dir = correctDir(edges->getDir(), edges->getComp());  // correctDir by MP mapping PE
				EdgeDir StartVertixID_dir;  // correctDir by MP mapping PE
				// Recursively add neighbors by distance
				addNeighborsToSubgraph_dfs(pRootVertex, pSubgraph, StartVertixID_dir, &stop, &Turn_off_flag, &final_distance, upper_distance, curr_distance, &leaf_count, &VertexCount, LargeDegreeVertexFromS, End_set, visit_Vertex_set);
				pSubgraph->writeDot("FINAL.dot");
			//	std::cout<<"dis="<<final_distance<<std::endl;
				if(LargeDegreeVertexFromS!="")
				{
					//std::cout<<"LDFS="<<LargeDegreeVertexFromS<<"\n";
					LDFS.push_back(LargeDegreeVertexFromS);
				}
			}
			if(stop==1)  // break of muti start
			{
				//std::cout << "mp" << mateNum+1 << "\t StartVertex: " << StartVertixID << " EndVertex: " << LargeDegreeVertexFromS << " upper_dis= " << upper_distance << std::endl;
				break;
			}
		}  // for start point
		
		if(Mode_Control==0 && final_distance==0)
			Mode_Control=1; // mode1: SE failed
		else if (final_distance!=0)
		{
			outStartVertexID=StartVertixID;
			//outEndVertexID=EndVertixID;
			Mode_Control=0;
		}
		VertexCount1=VertexCount;
		
		//for(size_t i =0 ;i<LDFS.size();i++)
			//std::cout<<"*****LDFS="<<LDFS[i]<<"\n";
			

		//
		//  Step.2 End to Start (EndVertex to LargeDegreeVertexFromS)
		if(Mode_Control==1)  // if Start to End Faild
		{
			//std::cout << "ESmod:::ESmode Step.2 " << std::endl;
			EdgeHashSet* End_set_ES = new EdgeHashSet;
			End_set_ES->set_empty_key("");

			//            if(LargeDegreeVertexFromS!="")
			//              End_set_ES->insert(LargeDegreeVertexFromS);    // End_point Hash_set
			if(LDFS.size()!=0)
			{
				for(size_t i=0;i<LDFS.size();i++)
				{
					if(LDFS[i]!="")
						End_set_ES->insert(LDFS[i]);
				}
			}
			//else if (LargeDegreeVertexFromS=="")
			else if (LDFS.size()==0)
			{
						std::cout << "mp" << mateNum+1 << "\t Step.2 LargeDegreeVertexFromS is NULL\n" ;
						continue;
			}

			EdgeHashSet* visit_Vertex_set_ES = new EdgeHashSet;
			visit_Vertex_set_ES->set_empty_key("");

			for(size_t i = 0; i < EndVertixID_list.size(); i++)    // for each start point // left to right search
			{
				if (i>10) break;	// only try three StartVertex
				StartVertixID=EndVertixID_list[i];
				pRootVertex = pGraph->getVertex(StartVertixID);	// Get the start(root) vertex from main graph
				curr_distance=0, final_distance=0;
				stop=0, Turn_off_flag=0;
				leaf_count=opt::LeafCount-50; 	// no necessary to DFS search , just find LargeDegreeVertexFromE
				VertexCount=0;
				LargeDegreeVertexFromE = "";	// no necessary initial because only keep one time
				upper_distance=opt::depth;

				if(pRootVertex == NULL)
				{
					std::cout << "mp" << mateNum+1 << "\t Step.2 StartVertex::" << StartVertixID << " not found in the graph.\n";
					continue;
				}
				else
				{
					EdgeDir StartVertixID_dir;  // correctDir by MP mapping PE
					//EdgeDir StartVertixID_dir = correctDir(edges_start[i]->getDir(), edges_start[i]->getComp());  // correctDir by MP mapping PE
					// Recursively add neighbors by distance
					addNeighborsToSubgraph_dfs(pRootVertex, pSubgraph, StartVertixID_dir, &stop, &Turn_off_flag, &final_distance, upper_distance, curr_distance, &leaf_count, &VertexCount, LargeDegreeVertexFromE, End_set_ES, visit_Vertex_set_ES);
					if(LargeDegreeVertexFromE!="")
						LDFE.push_back(LargeDegreeVertexFromE);
				}
				if(stop==1)  // break of muti start
				{
					//std::cout << "mp" << mateNum+1 << "StartVertex: " << StartVertixID << " EndVertex: " << LargeDegreeVertexFromS << " upper_dis= " << upper_distance << std::endl;
					break;
				}
			}  // for start point 

			delete End_set_ES;
			delete visit_Vertex_set_ES;

			VertexCount2=VertexCount;
		}
		if(Mode_Control==1 && final_distance!=0)
			Mode_Control=2; // mode2 = SE failed & ES sucessful
		else if(Mode_Control==1 && final_distance==0)
			Mode_Control=3; // mode3 = SE & ES failed 

		//
		//  Step.3 Start to End (StartVertex to LargeDegreeVertexFromS)
		//std::string outStartVertexID="",outEndVertexID="";
		if(Mode_Control==3)  // if Start to End Faild
		{
			//std::cout << "SEmod:::SEmode Step.3 " << std::endl;
			EdgeHashSet* End_set_SE = new EdgeHashSet;
			End_set_SE->set_empty_key("");
			// if(LargeDegreeVertexFromS!="")
			if(LDFS.size()!=0)
			{
				for(size_t i=0;i<LDFS.size();i++)
				{
					if(LDFS[i]!="")
					{
						//std::cout<<"LargeDegree>5= "<<LDFS[i]<<std::endl;
						End_set_SE->insert(LDFS[i]);  
					}
				}
}						// End_point Hash_set
			//else if (LargeDegreeVertexFromS=="")
			else if(LDFS.size()==0)
					{
						std::cout << "mp" << mateNum+1 << "\t Step.3 LargeDegreeVertexFromS is NULL\n" ;
						continue;
					}

			EdgeHashSet* visit_Vertex_set_SE = new EdgeHashSet;
			visit_Vertex_set_SE->set_empty_key("");
		//	std::cout<<"*********size_end="<<LDFS.size()<<"******Start="<<LDFS[0]<<std::endl;

			for(size_t i = 0; i < StartVertixID_list.size(); ++i)    // for each start point // Start to End search
				//for(size_t i = 0; i < 1; i++)
			{
				//std::cout<<"size= "<<StartVertixID_list.size()<<"\t StartVertixID_list[i]="<<StartVertixID_list[i]<<std::endl;
				if(i>10) break;
				StartVertixID=StartVertixID_list[i];
				Comp_flagS=0;
				pRootVertex = pGraph->getVertex(StartVertixID);	// Get the start(root) vertex from main graph
				curr_distance=0, final_distance=0;
				stop=0, Turn_off_flag=0;
				//leaf_count=opt::LeafCount;
				leaf_count=0;
				VertexCount=0;
				upper_distance=opt::depth;
				LargeDegreeVertexFromS="";

				if(pRootVertex == NULL)
				{
					std::cout << "mp" << mateNum+1 << "\t StartVertex::" << StartVertixID << " not found in the graph.\n";
					continue;
				}
				else
				{
					EdgeDir StartVertixID_dir;  // correctDir by MP mapping PE
					//EdgeDir StartVertixID_dir = correctDir(edges_start[i]->getDir(), edges_start[i]->getComp());  // correctDir by MP mapping PE
					// Recursively add neighbors by distance
		
					addNeighborsToSubgraph_dfs(pRootVertex, pSubgraph, StartVertixID_dir, &stop, &Turn_off_flag, &final_distance, upper_distance, curr_distance, &leaf_count, &VertexCount, LargeDegreeVertexFromS, End_set_SE, visit_Vertex_set_SE);
				
				}
				if(stop==1)  // break of muti start
				{
				//	std::cout << "Step.3 mp" << mateNum+1 << " Flag: "<< Comp_flagS << " StartVertex: " << StartVertixID << " EndVertex: " << LargeDegreeVertexFromS << std::endl;
					outStartVertexID=StartVertixID;
					break;
				}
			}  // for start point 

			delete End_set_SE;
			delete visit_Vertex_set_SE;

			VertexCount3=VertexCount;
		}

		//
		//  Step.4 End to Start (EndVertex to LargeDegreeVertexFromL)
		//
		if(Mode_Control==3)  // if Start to End Faild
		{
			//std::cout << "ESmod:::ESmode Step.4 " << std::endl;
			EdgeHashSet* End_set_ES = new EdgeHashSet;
			End_set_ES->set_empty_key("");
			//if(LargeDegreeVertexFromE!="")
			if(LDFE.size()!=0)
			{
				for(size_t i=0;i<LDFE.size();i++)
				{
					if(LDFE[i]!="")
						End_set_ES->insert(LDFE[i]);  
				}
			}						// End_point Hash_set
			// else if (LargeDegreeVertexFromE=="")
			else if(LDFE.size()==0)
			{
				std::cout << "mp" << mateNum+1 << "\t Step.4 LargeDegreeVertexFromE is NULL\n" ;
				continue;
			}

			EdgeHashSet* visit_Vertex_set_ES = new EdgeHashSet;
			visit_Vertex_set_ES->set_empty_key("");
//std::cout<<"*********size_end="<<LDFE.size()<<"******End="<<LDFE[0]<<std::endl;
			for(size_t i = 0; i < EndVertixID_list.size(); ++i)    // for each start point // left to right search
				//for(size_t i = 0; i < 1; i++)    // for each start point // left to right search
			{
				//std::cout<<"fuck you\n";
				if(i>10) break;
				StartVertixID=EndVertixID_list[i];  // Get the start(root) vertex from main graph
				//std::cout<<"StartVertixID="<<StartVertixID<<std::endl;
				Comp_flagE=0;
				pRootVertex = pGraph->getVertex(StartVertixID);
				curr_distance=0, final_distance=0;
				stop=0, Turn_off_flag=0;
				//leaf_count=opt::LeafCount;
				leaf_count=0;
				VertexCount=0;
				upper_distance=opt::depth;

				if(pRootVertex == NULL)
				{
					std::cout << "mp" << mateNum+1 << "\t Step.4 StartVertex " << StartVertixID << " not found in the graph.\n";
					continue;
				}
				else
				{
					EdgeDir StartVertixID_dir ;  // correctDir by MP mapping PE
					//EdgeDir StartVertixID_dir = correctDir(edges_start[i]->getDir(), edges_start[i]->getComp());  // correctDir by MP mapping PE
					// Recursively add neighbors by distance
					addNeighborsToSubgraph_dfs(pRootVertex, pSubgraph, StartVertixID_dir, &stop, &Turn_off_flag, &final_distance, upper_distance, curr_distance, &leaf_count, &VertexCount, LargeDegreeVertexFromE, End_set_ES, visit_Vertex_set_ES);
					//std::cout<<"StartVertixID="<<StartVertixID<<std::endl;
				}
				if(stop==1)  // break of muti start
				{
					//std::cout << "Step.4 mp" << mateNum+1 << " Flag: "<< Comp_flagE << " StartVertex: " << StartVertixID << " EndVertex: " << LargeDegreeVertexFromE << std::endl;
					outEndVertexID=StartVertixID;
					break;
				}
			}  // for start point 

			delete End_set_ES;
			delete visit_Vertex_set_ES;

			VertexCount4=VertexCount;
		}

		//
		// OutPut: simplify contigs_graph and output contigs
		//int thread_num=omp_get_thread_num();
		//std::cout << "mp = " <<mateNum+1 <<"\tfuck"<< std::endl;
		//std::cout <<"Mode="<<Mode_Control<< "\t mp=" <<mateNum+1<< "\t outStartVertexID= " << outStartVertexID<<"\t outEndVertexID"<<outEndVertexID<< std::endl;
		//pSubgraph->writeDot("Over.dot");
		//OutputContig(pSubgraph, Mode_Control, mateNum, final_distance, thread_num, outStartVertexID, outEndVertexID, Comp_flagS, Comp_flagE);
		

		delete visit_Vertex_set;
		delete End_set;
		delete pSubgraph;
		

		m_end = clock();
		std::cout << "\tTotal Time:" << ((m_end-m_start)/(double)(CLOCKS_PER_SEC)) ; //<< std::endl; 
		std::cout << " 1:" << VertexCount1 << " 2:" << VertexCount2 << " 3:" << VertexCount3 << " 4:" << VertexCount4 << std::endl;
*/
	} // for mateNum (each mate pair)

	delete pGraph;
	delete pBWT;
	delete pRBWT;
	delete pSSA;
	delete pSSA2;
	delete pIntervalCache;
	delete pTargetReadInfoTable;

}
///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
void addNeighborsToSubgraph_dfs(Vertex* pCurrVertex, StringGraph* pSubgraph, EdgeDir dir, int *stop, int *Turn_off_flag, unsigned int *final_distance, unsigned int upper_distance, unsigned int curr_distance, unsigned int *leaf_count, unsigned int *VertexCount, std::string &LargeDegreeVertex, EdgeHashSet* End_set, EdgeHashSet* visit_Vertex_set)
{

	// These are the edges in the main graph
	EdgePtrVec edges=pCurrVertex->getEdges(dir); 
//std::cout << " vertex:" << pCurrVertex->getID()  << " curr_dis:" << curr_distance << " edge_num:" << edges.size() << std::endl;
//std::cout<<"******stop="<<*stop<<std::endl;
	*VertexCount=*VertexCount+1;

	//==================================================================================
	// search control secion
	//==================================================================================

	if ( LargeDegreeVertex == "" && edges.size() >= 3 && *VertexCount>1)   // Keep LargeDegreeVertex
	{
		LargeDegreeVertex = pCurrVertex->getID();
		//std::cout << "LargeDegreeVertex: " << LargeDegreeVertex << std::endl;
	}

	if(opt::TotalVertexCount > 0 && *VertexCount > opt::TotalVertexCount)  // search control (Vertex number)
	{
		//std::cout << "TouchVertexCount: " << *VertexCount << std::endl;
		*Turn_off_flag=1;
		return;
	}	

	if(curr_distance > upper_distance)  // search control (depth)
	{
		*leaf_count=*leaf_count+1;
		//std::cout<<"curr="<<curr_distance<<std::endl;
		//std::cout << "leaf_count(depth): " << *leaf_count << std::endl;
		return;
	}	

	if(opt::LeafCount > 0 && *leaf_count > opt::LeafCount) // search control (leaf)
	{
		//std::cout << "out of leaf_count(up_bound): " << *leaf_count << std::endl;
		*Turn_off_flag=1;
		return;
	}

	if(edges.size()==0)  // if vertex degree=0 leaf_count++
	{
		*leaf_count=*leaf_count+1;
		//std::cout << "leaf_count(degree): " << *leaf_count << std::endl;
		return;
	}

	//==================================================================================
	//==================================================================================

	// DFS start 
	// DFS start 
	for(size_t i = 0; i < edges.size(); ++i)
	{
		Vertex* pY = edges[i]->getEnd();
		EdgeHashSetIter EHSiter = visit_Vertex_set->find(pY->getID());

		if( !(EHSiter != visit_Vertex_set->end()) && *stop!=1 )
		{ 
			visit_Vertex_set->insert(pY->getID());  // insert visited vertex to edge_hash
			// Calculated distance
			curr_distance = curr_distance + (edges[i])->getEnd()->getSeqLen() - (edges[i])->getOverlap().getOverlapLength(0);

			//printf area
			//printf area
			//std::cout << " vertex:" << pCurrVertex->getID()  << " curr_dis:" << pCurrVertex->getDistance() << " edge_num:" << edges.size() << std::endl;
			/*      
				if(edges.size()>1)
				{
				for (unsigned int ijk=0 ; ijk<edges.size() ; ijk++)
				{
				if(ijk<3)
				std::cout << "          " << (edges[ijk])->getOverlap().getOverlapLength(0) << " VertixID: " << (edges[ijk])->getEnd()->getID() << std::endl;
				}
				}
			 */    
			//printf area
			//printf area

			// End point hash 
			EdgeHashSetIter RHSiter = End_set->find((edges[i])->getEnd()->getID()); // find next vertex in EndSet_hash or not.

			if (RHSiter != End_set->end())   // find EndVertex successful
			{
				//std::cout<<edges[i]->getEnd()->getID()<<std::endl;
				*stop=1;
				//*Turn_off_flag=1;
				*final_distance=curr_distance;
				//std::cout<<"stop="<<*stop<<std::endl;
				//std::cout << " vertex:" << pCurrVertex->getID()  << " curr_dis:" << curr_distance << " edge_num:" << edges.size() << std::endl;
				copyVertexToSubgraph(pSubgraph, pCurrVertex);
				copyVertexToSubgraph(pSubgraph, pY);
				Overlap ovr = edges[i]->getOverlap();
				SGAlgorithms::createEdgesFromOverlap(pSubgraph, ovr, true);
				//std::cout<<"************* \n";
				return;
			}
			if(*stop==1)
		{
			//std::cout<<"stop="<<*stop<<std::endl;
			Vertex* pZ = edges[i]->getEnd();
			//std::cout << " vertex:" << pCurrVertex->getID()  << " curr_dis:" << curr_distance << " edge_num:" << edges.size() << std::endl;
			copyVertexToSubgraph(pSubgraph, pCurrVertex);
			copyVertexToSubgraph(pSubgraph, pZ);
			Overlap ovr = edges[i]->getOverlap();
			SGAlgorithms::createEdgesFromOverlap(pSubgraph, ovr, true);
			return;
		}
			else if(*stop!=1)
			{
			// Correct dir for complement edge
			assert(edges[i]->getDir() == dir);
			EdgeDir corrected_dir = correctDir(edges[i]->getDir(), edges[i]->getComp()); 
			addNeighborsToSubgraph_dfs(pY, pSubgraph, corrected_dir, stop, Turn_off_flag, final_distance, upper_distance, curr_distance, leaf_count, VertexCount, LargeDegreeVertex, End_set, visit_Vertex_set);
			}
			if(*Turn_off_flag==1)	// set by leaf_count or vertex_count
				break;
		}
		/*
		   else if (EHSiter != visit_Vertex_set->end()) // if the vertex has been travaled leaf_count++;
		   {
		 *leaf_count=*leaf_count+1;	
		//std::cout << "leaf_count(visited): " << *leaf_count << std::endl;
		}
		 */		

		// Reached EndPoint :: Create subgraph to output Reads
		/*if(*stop==1)
		{
			//std::cout<<"stop="<<*stop<<std::endl;
			Vertex* pY = edges[i]->getEnd();
			//std::cout << " vertex:" << pCurrVertex->getID()  << " curr_dis:" << curr_distance << " edge_num:" << edges.size() << std::endl;
			copyVertexToSubgraph(pSubgraph, pCurrVertex);
			copyVertexToSubgraph(pSubgraph, pY);
			Overlap ovr = edges[i]->getOverlap();
			SGAlgorithms::createEdgesFromOverlap(pSubgraph, ovr, true);
			return;
		}*/
		//Reached EndPoint :: Create subgraph to output Reads
	}
}

void addNeighborsToSubgraph(Vertex* pCurrVertex,StringGraph* pSubgraph,int span)
{
if(span <=0)
	return ;
EdgePtrVec edges = pCurrVertex->getEdges();	



}
void OutputContig_bfs_case1(StringGraph* pSubgraph,unsigned int mateNum,unsigned int total_distance,int thread_num,std::string outStartVertexID)//for case 1
{
SGTrimVisitor trimVisit(opt::trimLengthThreshold);
std::string string_thread_num;
string_thread_num = int2str(thread_num);
std::string outInnerReads = "mp" + string_thread_num + "-InnerReads.fa";  // tt = threads_num
//std::string outOuterReads = "mp" + string_thread_num + "-OuterReads.fa";
//std::string outMode2Reads = "mp" + string_thread_num + "-EndReads.fa";
//std::string outMode3Reads = "mp" + string_thread_num + "-SeparateReads.fa";
char char_mpNum[20];
	sprintf(char_mpNum,"%d",mateNum+1);
	std::string string_mpNum = char_mpNum ; 
	std::string ContigName = "mp" + string_mpNum + "_";
	if(opt::numTrimRounds>0)
	{
		int numTrims = opt::numTrimRounds;
		while(numTrims-- > 0)
			pSubgraph->visit(trimVisit);
	}
//12-22	pSubgraph->simplify_MP(outStartVertexID);
	pSubgraph->simplify();
	//assert( pSubgraph->getAllVertices().size() == 1);
        std::cout << pSubgraph->getFirstVertex()->getStr();	
        ContigName = "mp" + string_mpNum + "_case1_";
	SGMPFastaVisitor av(outInnerReads);
	pSubgraph->visit(av);
	std::cout << "mp" << mateNum+1 << "_case1" << "\t Find a path. final_dis= " << total_distance<<std::endl ;

}/*
void OutputContig_bfs_case2(StringGraph* pSubgraph_start,StringGraph* pSubgraph_end,unsigned int mateNum,int thread_num,std::string outStartVertexID,std::string outEndVertexID)//for case 2
{
SGTrimVisitor trimVisit(opt::trimLengthThreshold);
std::string string_thread_num;
string_thread_num = int2str(thread_num);
//std::string outInnerReads = "mp" + string_thread_num + "-InnerReads.fa";  // tt = threads_num
//std::string outOuterReads = "mp" + string_thread_num + "-OuterReads.fa";
//std::string outMode2Reads = "mp" + string_thread_num + "-EndReads.fa";
std::string outMode3Reads = "mp" + string_thread_num + "-SeparateReads.fa";
char char_mpNum[20];
	sprintf(char_mpNum,"%d",mateNum+1);
	std::string string_mpNum = char_mpNum ; 
	std::string ContigName = "mp" + string_mpNum + "_";
	if(opt::numTrimRounds>0)
	{
		int numTrims = opt::numTrimRounds;
		while(numTrims-- > 0)
		{
			pSubgraph_start->visit(trimVisit);
			pSubgraph_end->visit(trimVisit);
		}	
	}
	pSubgraph_start->simplify_MP(outStartVertexID);
	pSubgraph_end->simplify_MP(outEndVertexID);
	ContigName = "mp" + string_mpNum + "_case2_";
	SGMPFastaVisitor av(outMode3Reads);
	pSubgraph_start->visit(av);
	pSubgraph_end->visit(av);
	std::cout << "mp" << mateNum+1 << "_case2 \n";

}*/
//
//	output 
//
/*
void OutputContig(StringGraph* pSubgraph, int Mode_Control, unsigned int mateNum, unsigned int final_distance, int thread_num, std::string outStartVertexID, std::string outEndVertexID, int Comp_flagS, int Comp_flagE)
{
	//std::cout<<"mode="<<Mode_Control<<"\t dis="<<final_distance<<std::endl;
	// Simplify and rename SuperReads
	//std::cout << "Output mp=" <<mateNum+1<< "\t outStartVertexID= " << outStartVertexID<<"\t outEndVertexID"<<outEndVertexID<< std::endl;
	SGTrimVisitor trimVisit(opt::trimLengthThreshold);

	//Initial outputFileName
	std::string string_thread_num;
	string_thread_num = int2str(thread_num);
	std::string outInnerReads = "mp" + string_thread_num + "-InnerReads.fa";  // tt = threads_num
	std::string outOuterReads = "mp" + string_thread_num + "-OuterReads.fa";
	std::string outMode2Reads = "mp" + string_thread_num + "-EndReads.fa";
	std::string outMode3Reads = "mp" + string_thread_num + "-SeparateReads.fa";

	char char_mpNum[20];
	sprintf(char_mpNum,"%d",mateNum+1);
	std::string string_mpNum = char_mpNum ; 
	std::string ContigName = "mp" + string_mpNum + "_";
//std::cout<<"mode="<<Mode_Control<<"\t dis="<<final_distance<<std::endl;
	// MatePair insert size variance control
	// innerReads ( min_distance <= superRead length <= max_distance)
	// outerReads ( superRead length > max_distance OR superRead length < min_distance)
	unsigned int max_distance = opt::mean_d + opt::var_d ;  // int max_distance = mean + variance ;
	unsigned int min_distance = opt::mean_d - opt::var_d ;  // int min_distance = mean - variance ;
//std::cout<<"Max="<<max_distance<<std::endl;
	//std::cout<<"mode="<<Mode_Control<<"\t dis="<<final_distance<<std::endl;
	//pSubgraph->writeDot("test.dot");
	if(Mode_Control==0 && final_distance<=max_distance && final_distance>=min_distance)
	{	
		// Remove dead-end branches from the graph
		if(opt::numTrimRounds > 0) 
		{
			int numTrims = opt::numTrimRounds;
			while(numTrims-- > 0)
				pSubgraph->visit(trimVisit);
		}
//std::cout<<"Start="<<outStartVertexID<<"\t End="<<outEndVertexID<<std::endl;
		
		pSubgraph->simplify_MP(outStartVertexID, outEndVertexID);

		int ContigNum = pSubgraph->simplify_CountVertex();
		//std::cout<<"Num="<<ContigNum<<std::endl;
		if(ContigNum==1)
		{
			ContigName = "mp" + string_mpNum + "_mode0_" ;
			pSubgraph->renameVertices_MP(ContigName, outStartVertexID, outEndVertexID, Comp_flagS, Comp_flagE);
			SGMPFastaVisitor av(outInnerReads);
			pSubgraph->visit(av);
			std::cout << "mp" << mateNum+1 << "_mode0" << "\t Find a path. final_dis= " << final_distance ;
		}
	}
	else if (Mode_Control==0 && final_distance!= 0 && ( final_distance > max_distance || final_distance < min_distance ) )
	{
		//std::cout<<"Start="<<outStartVertexID<<"\t End="<<outEndVertexID<<std::endl;
		pSubgraph->simplify_MP(outStartVertexID, outEndVertexID);

		int ContigNum = pSubgraph->simplify_CountVertex();
		if(ContigNum==1)
		{
			ContigName = "mp" + string_mpNum + "_mode0_" ;
			//std::cout<<"start="<<outStartVertexID<<"\t end="<<outEndVertexID<<std::endl;
			pSubgraph->renameVertices_MP(ContigName, outStartVertexID, outEndVertexID, Comp_flagS, Comp_flagE);
			SGMPFastaVisitor av(outOuterReads);
			pSubgraph->visit(av);
			std::cout << "mp" << mateNum+1 << "_mode0" << "\t Out of dis. real_dis=" << final_distance ;
		}
		else
		{
			std::cout << "mp" << mateNum+1 << "_mode0" << "\t Contig Num. Strange.";
			//pSubgraph->writeDot("test.dot");
		}
	}
	else if ( Mode_Control==2 && final_distance!= 0 )
	{
		pSubgraph->simplify_MP(outStartVertexID, outEndVertexID);	
		int ContigNum = pSubgraph->simplify_CountVertex();
		if(ContigNum==1)
		{
			ContigName = "mp" + string_mpNum + "_mode2_" ; 
			pSubgraph->renameVertices_MP(ContigName, outStartVertexID, outEndVertexID, Comp_flagS, Comp_flagE);
			SGMPFastaVisitor av(outMode2Reads);
			pSubgraph->visit(av);
			std::cout << "mp" << mateNum+1 << "_mode2" ;
		}
		else
			std::cout << "mp" << mateNum+1 << "_mode2" << "\t Contig Num. Strange.";
	}
	else if ( Mode_Control==3 && final_distance!= 0 )
	{
		//std::cout<<"mode="<<Mode_Control<<"\t dis="<<final_distance<<std::endl;
		//std::cout<<"mode="<<Mode_Control<<std::endl;
		std::string FlagS,FlagE;
		//std::cout<<"mode="<<Mode_Control<<"\t dis="<<final_distance<<std::endl;
		FlagS=int2str(Comp_flagS);
		//std::cout<<"mode="<<Mode_Control<<"\t dis="<<final_distance<<std::endl;
		FlagE=int2str(Comp_flagE);
		//std::cout<<"mode="<<Mode_Control<<"\t dis="<<final_distance<<std::endl;
		//std::cout << "Output mp=" <<mateNum+1<< "\t outStartVertexID= " << outStartVertexID<<"\t outEndVertexID"<<outEndVertexID<< std::endl;
pSubgraph->simplify_MP(outStartVertexID, outEndVertexID);//problem????
//pSubgraph->simplify();
		//std::cout<<"mode="<<Mode_Control<<"\t dis="<<final_distance<<std::endl;
//std::cout<<"mode="<<Mode_Control<<std::endl;
		int ContigNum = pSubgraph->simplify_CountVertex();
		if(ContigNum==2)
		{
			ContigName = "mp" + string_mpNum + "_mode3_" + FlagS + ":" + FlagE + "_" ;
			pSubgraph->renameVertices_MP(ContigName, outStartVertexID, outEndVertexID, Comp_flagS, Comp_flagE);
			SGMPFastaVisitor av(outMode3Reads);
			pSubgraph->visit(av);
			std::cout << "mp" << mateNum+1 << "_mode3" ;
		}
		else
			std::cout << "mp" << mateNum+1 << "_mode3" << "\t Contig Num. Strange.";
	}
	else
	{
		std::cout << "mp" << mateNum+1 << "\t Can't find a path " <<"\tMode_Control="<<Mode_Control<<"\tdis="<<final_distance<<"\n";
	}
}
*/
// 
// Handle command line arguments
//
void parseSubgraphOptions(int argc, char** argv)
{
	// Set defaults
	std::string prefix = "default" ;
	bool die = false;
	for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) 
	{
		std::istringstream arg(optarg != NULL ? optarg : "");
		switch (c) 
		{
			//case 'i': opt::bIsInptList = true; break;
			case 'o': arg >> prefix; break;
			case 'd': arg >> opt::Edir; break;
			case 's': arg >> opt::span; break;
			case 'm': arg >> opt::mean_d; break;
			case 'a': arg >> opt::var_d; break;
			case 't': arg >> opt::depth; break;
			case 'F': arg >> opt::LeafCount; break;
			case 'V': arg >> opt::TotalVertexCount; break;
			case 'L': arg >> opt::trimLengthThreshold; break;
			case 'X': arg >> opt::numTrimRounds; break;
			case '?': die = true; break;
			case 'v': opt::verbose++; break;
			case OPT_TR: opt::bPerformTR = true; break;
			case OPT_FQ: opt::loadFastQ = true; break;
			case OPT_FA: opt::loadFastA = true; break;

			case OPT_HELP:
				     std::cout << SUBGRAPH_USAGE_MESSAGE;
				     exit(EXIT_SUCCESS);
			case OPT_VERSION:
				     std::cout << SUBGRAPH_VERSION_MESSAGE;
				     exit(EXIT_SUCCESS); 
		}
	}

	opt::outSubgraph      = prefix + "-subgraph.asqg";    // opt::outSubgraph = "subgraph.asqg.gz";
	opt::outContigsFile   = prefix + "-contig.fa";
	//opt::outGraphFile     = prefix + "-contig-graph.asqg";
	//opt::outSimplifyGraph = prefix + "rmTranCont-graph.asqg";

	if (argc - optind <	3) 
	{
		std::cerr << SUBPROGRAM ": missing arguments\n";
		die = true;
	} 
	else if (argc - optind > 3) 
	{
		std::cerr << SUBPROGRAM ": too many arguments\n";
		die = true;
	}

	if (die) 
	{
		std::cout << "\n" << SUBGRAPH_USAGE_MESSAGE;
		exit(EXIT_FAILURE);
	}

	// Parse the index file
	opt::indexPrefix = argv[optind++];
	// Parse the mate-pair file
	opt::MpFile = argv[optind++];

	if(opt::indexPrefix.empty())
	{
		std::cerr << SUBPROGRAM ": missing indexPrefix\n";
		exit(EXIT_FAILURE);
	}
	if(opt::MpFile.empty())
	{
		std::cerr << SUBPROGRAM ": missing MpFile\n";
		exit(EXIT_FAILURE);
	}
	// Parse the input filename
	opt::asqgFile  = argv[optind++];
	if(opt::asqgFile.empty())
	{
		std::cerr << SUBPROGRAM ": missing ASQGFile\n";
		exit(EXIT_FAILURE);
	}

} 
//
// Load *.fastq to locate Start and End
//
void load_MpSeqInfo()
{
	//	load *.fastq
	int count_line=0;
	std::string str;
	//std::string file1 = "sga_pre_mp.fastq";	
	std::string file1 = opt::MpFile;
	std::ifstream fin(file1.c_str(),std::ios::in);
	if (fin == NULL)
	{
		std::cout << " Can't open " << file1 << std::endl;
		exit(0);
	}

	if(opt::loadFastQ)
	{
		while (getline(fin, str))
		{
			if(count_line%8==0)
				opt::MpInfo1_ID.push_back(str);
			else if(count_line%8==1)
				opt::MpInfo1_SEQ.push_back(str);
			else if (count_line%8==4)
				opt::MpInfo2_ID.push_back(str);
			else if(count_line%8==5)
				opt::MpInfo2_SEQ.push_back(str);

			count_line++;
		}
		fin.close();
	}
	else if(opt::loadFastA)
	{
		while (getline(fin, str))
		{
			if(count_line%5==0)
				opt::MpInfo1_ID.push_back(str);
			else if(count_line%5==1)
				opt::MpInfo1_SEQ.push_back(str);
			else if (count_line%5==2)
				opt::MpInfo2_ID.push_back(str);
			else if(count_line%5==3)
				opt::MpInfo2_SEQ.push_back(str);
                        else if(count_line%5==4)
                                opt::MpInfo_length.push_back(str);
			count_line++;
		}
		fin.close();
	}
	else 
	{
		std::cout << " Can't get maFile format (--fastq/--fasta) " << file1 << std::endl;
		exit(0);
	}
}
/*
   void load_MpSeqInfo_fasta()
   {
//	load *.fastq
int count_line=0,fq_flag=0;
std::string str;
//std::string file1 = "sga_pre_mp.fastq";	
std::string file1 = opt::MpFile;
std::ifstream fin(file1.c_str(),std::ios::in);
if (fin == NULL)
{
std::cout << " Can't open " << file1 << std::endl;
exit(0);
}

if(fq_flag==0)
{
while (getline(fin, str))
{
if(count_line%4==0)
opt::MpInfo1_ID.push_back(str);
else if(count_line%4==1)
opt::MpInfo1_SEQ.push_back(str);
else if (count_line%4==2)
opt::MpInfo2_ID.push_back(str);
else if(count_line%4==3)
opt::MpInfo2_SEQ.push_back(str);

count_line++;
}
fin.close();
}
}
 */
// We do not directly add the vertex pointer to the graph, as each graph
// manages its set of vertices and we do not want to double-free the vertices
void copyVertexToSubgraph(StringGraph* pSubgraph, const Vertex* pVertex)
{
	// Make sure the vertex hasn't been added yet
	//std::cout<<"shit"<<std::endl;
	if(pSubgraph->getVertex(pVertex->getID()) == NULL)
	{
		Vertex* pCopy = new(pSubgraph->getVertexAllocator()) Vertex(pVertex->getID(), pVertex->getSeq().toString());
		pSubgraph->addVertex(pCopy);
	}
}

// int to string
std::string int2str(int &i)   // int to tring
{
	std::string s;
	std::stringstream ss(s);
	ss << i;
	return ss.str();
}



