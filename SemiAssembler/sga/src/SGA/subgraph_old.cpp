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
void addNeighborsToSubgraph_dfs(Vertex* pCurrVertex, StringGraph* pSubgraph, EdgeDir dir, int *stop, unsigned int *final_distance, unsigned int upper_distance, unsigned int curr_distance, unsigned int *leaf_count, std::string &LargeDegreeVertex, EdgeHashSet* End_set, EdgeHashSet* visit_Vertex_set, unsigned int *VertexCount);
void load_MpSeqInfo();  // load *.fastq
void OutputContig(StringGraph* pSubgraph, int Mode_Control, unsigned int mateNum, unsigned int final_distance, int thread_num, std::string outStartVertixID, std::string outEndVertixID, int Comp_flagS, int Comp_flagE);
void copyVertexToSubgraph(StringGraph* pSubgraph, const Vertex* pVertex);

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
	static unsigned int mean_d = 4500;         // mean of insert size
	static unsigned int var_d = 1500;          // variance of insert size

	// Search Control parameters
	static int depth = mean_d + var_d + 200;   // subgraph depth
	static unsigned int LeafCount = 0;
	static unsigned int VertexCount = 0;

	// Trim parameters
	static int numTrimRounds = 5;
	static unsigned int trimLengthThreshold = 300;

	// Vistor Control
	static bool bPerformTR  = false;            // SGTransitiveReductionVisitor

	// Data array
	static std::vector<std::string> MpInfo1_ID;     // mp_ID1_array
	static std::vector<std::string> MpInfo1_SEQ;    // mp_ID1_array
	static std::vector<std::string> MpInfo2_ID;     // mp_ID2_array
	static std::vector<std::string> MpInfo2_SEQ;    // mp_ID2_array


}

static const char* shortopts = "o:s:d:m:a:t:F:L:X:V:";

enum { OPT_HELP = 1, OPT_VERSION, OPT_TR };

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
	//#pragma omp parallel for num_threads(60) private(mateNum) 
	//===========================================================
	for(mateNum=0 ; mateNum < mateNum_up_bound ; mateNum++)  // for each mate pair
	{	
		clock_t tt_start, tt_end;
		tt_start = clock();

		clock_t m_start,m_end;	//Timer* pTimerE = new Timer("Each matepair"); 
		m_start = clock();

		std::string StartVertixID = opt::MpInfo1_ID[mateNum];  // start ID in mp_ID1 (xxx/1)
		std::string EndVertixID   = opt::MpInfo2_ID[mateNum];  // end   ID in mp_ID2 (xxx/2)

		std::vector<std::string> StartVertixID_list;
		std::vector<std::string> EndVertixID_list;

		//
		// Get EndVertexID and insert to EndVetex_hash
		int substr_len = opt::MpInfo2_SEQ[mateNum].length();
		std::string kmer = opt::MpInfo2_SEQ[mateNum].substr(0,substr_len-1);   // EndVertex_seq
BWTInterval interval;
                for(int i=0;i<70;i++)
{
 kmer = opt::MpInfo2_SEQ[mateNum].substr(0,substr_len-1-i);
 interval = BWTAlgorithms::findInterval(indexSet, kmer);
if(interval.size()>0)
 break;
}


		//BWTInterval interval = BWTAlgorithms::findInterval(indexSet, kmer);	 // Find FM index interval for EndVertex_SEQ
		kmer = reverseComplement(kmer);
		BWTInterval interval_rev = BWTAlgorithms::findInterval(indexSet, kmer);	// Find FM index interval for EndVertex_SEQ_reverse

		EdgeHashSet* End_set = new EdgeHashSet;  // EdgeHash for EndVertex (Start to End)
		End_set->set_empty_key("");

		//if( ( interval.size() + interval_rev.size() ) > 0 )		// MP hit num != 0
		if( interval.size() > 0 )		// MP hit num != 0
		{
			int count=0;
			for (int64_t j=interval.lower ; j<=interval.upper ; ++j)
			{
				if(count>5) break; // only keep 5 readID
				int64_t saIdx = j;
				const ReadInfo& targetInfo = pTargetReadInfoTable->getReadInfo(pSSA2->calcSA(saIdx,pBWT).getID()); // get PE id
				EndVertixID=targetInfo.id;
				End_set->insert(EndVertixID.c_str());    // insert to End_point Hash_set
				EndVertixID_list.push_back(EndVertixID); // push to EndVertex_list
				count++;
			}
			/*	for (int64_t j=interval_rev.lower ; j<=interval_rev.upper ; ++j)
				{
				int64_t saIdx = j;
				const ReadInfo& targetInfo = pTargetReadInfoTable->getReadInfo(pSSA2->calcSA(saIdx,pBWT).getID());
				EndVertixID=targetInfo.id;
				End_set->insert(EndVertixID.c_str());    // End_point Hash_set
				EndVertixID_list.push_back(EndVertixID);
				} */
		}
		else 
		{
			std::cout << "mp" << mateNum+1 << "\t End Vertex not found.\n";
			continue;
		}

		//
		// Get StartVertexID list 
		substr_len = opt::MpInfo1_SEQ[mateNum].length();
		//kmer = opt::MpInfo1_SEQ[mateNum].substr(0,substr_len-1);   // StartVertex_seq
//by little8

                for(int i=0;i<70;i++)
{
 kmer = opt::MpInfo1_SEQ[mateNum].substr(0,substr_len-1-i);
interval = BWTAlgorithms::findInterval(indexSet, kmer);
if(interval.size()>0)
 break;
}
//by little8

//����by little8
		//interval = BWTAlgorithms::findInterval(indexSet, kmer);	 // Find FM index interval for StartVertex_SEQ (kmer)
		kmer = reverseComplement(kmer);
		interval_rev = BWTAlgorithms::findInterval(indexSet, kmer);	// Find FM index interval for StartVertex_SEQ_reverse (kmer_reverse)

		//if( ( interval.size() + interval_rev.size() ) > 0 )		// MP hit num != 0
		if( interval.size() > 0 )		// MP hit num != 0
		{
			int count=0;
			for (int64_t j=interval.lower ; j<=interval.upper ; ++j)
			{
				if(count>5) break;	// only keep 5 readID
				int64_t saIdx = j;
				const ReadInfo& targetInfo = pTargetReadInfoTable->getReadInfo(pSSA2->calcSA(saIdx,pBWT).getID()); // get PE id
				StartVertixID=targetInfo.id;
				StartVertixID_list.push_back(StartVertixID);	// push to StartVertex_list
				count++;
			}
			/*	for (int64_t j=interval_rev.lower ; j<=interval_rev.upper ; ++j)
				{
				int64_t saIdx = j;
				const ReadInfo& targetInfo = pTargetReadInfoTable->getReadInfo(pSSA2->calcSA(saIdx,pBWT).getID());
				StartVertixID=targetInfo.id;
				StartVertixID_list.push_back(StartVertixID);
				}*/
		}
		else 
		{
			std::cout << "mp" << mateNum+1 << "\t Start Vertex not found.\n";
			continue;
		}

		tt_end = clock();
		//std::cout << "=====Time1:Start & End VertexList:" << ((tt_end-tt_start)/(double)(CLOCKS_PER_SEC)) << std::endl;
		tt_start = clock();


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
		unsigned int final_distance=0,upper_distance=0, curr_distance=0;
		unsigned int leaf_count=0;
		std::string LargeDegreeVertexFromS="";
		std::string LargeDegreeVertexFromE="";
		int Mode_Control=0;   // mode0: SE successful, mode2 = SE failed but ES sucessful, mode3 = SE & ES failed (convert to Nver SuperReads)
		int Comp_flagS=-1,Comp_flagE=-1;     // Reads Direction Check

		unsigned int VertexCount=0;
		unsigned int VertexCount1=0;
		unsigned int VertexCount2=0;
		unsigned int VertexCount3=0;
		unsigned int VertexCount4=0;

		//
		//	Start search matepair path
		//	Step1. StartVertex to EndVertex (if failed, goto Step2.)					// S--->Vs------------>E	Vs=LargeDegreeVertexFromS  (keep the Vertex that degree>5)
		//	Step2. EndVertex to LargeDegreeVertexFromS (if failed, goto Step3 & 4)		//		Vs<------Ve<---E	Ve=LargeDegreeVertexFromE
		//	Step3. StartVertex to LargeDegreeVertexFromS								// S--->Vs				
		//	Step4. EndVertex to LargeDegreeVertexFromE									//		         Ve<---E	
		//	fill the gap by N (by other program)										// S--->Vs  NNN  Ve<---E	fill the gap by N (by other program)
		//
		//  Step1. Start to End (Left to Right)
		for(size_t i = 0; i < StartVertixID_list.size(); i++)    // for each start point // Start to End search
		{
			//std::cout << "SEmod:::SEmode Step.1 " << std::endl;
			if(i>3) break;
			StartVertixID=StartVertixID_list[i];	// Get the start(root) vertex from main graph
			pRootVertex = pGraph->getVertex(StartVertixID);
			curr_distance=0;
			stop=0, final_distance=0, leaf_count=0;
			//LargeDegreeVertexFromS = "";	// no necessary initial because only keep one time
			upper_distance=opt::depth;

			if(pRootVertex == NULL)
			{
				std::cout << "mp" << mateNum+1 << "\t Step.1 StartVertex:: " << StartVertixID << " not found in the graph.\n";
				continue;
			}
			else
			{
				//EdgeDir StartVertixID_dir = correctDir(edges_start[i]->getDir(), edges_start[i]->getComp());  // correctDir by MP mapping PE
				EdgeDir StartVertixID_dir;  // correctDir by MP mapping PE
				// Recursively add neighbors by distance
				addNeighborsToSubgraph_dfs(pRootVertex, pSubgraph, StartVertixID_dir, &stop, &final_distance, upper_distance, curr_distance, &leaf_count, LargeDegreeVertexFromS, End_set, visit_Vertex_set, &VertexCount);
			}
			if(stop==1)  // break of muti start
			{
				//std::cout << "mp" << mateNum+1 << "StartVertex: " << StartVertixID << " EndVertex: " << LargeDegreeVertexFromS << " upper_dis= " << upper_distance << std::endl;
				break;
			}
		}  // for start point
		if(Mode_Control==0 && final_distance==0)
			Mode_Control=1; // mode1: SE failed
		else if (final_distance!=0)
			Mode_Control=0;

		VertexCount1=VertexCount;
		//std::cout << "   Step1.VertexCount=" << VertexCount << std::endl;
		tt_end = clock();
		//std::cout << "=====Time2:Step1:" << ((tt_end-tt_start)/(double)(CLOCKS_PER_SEC)) << std::endl;
		tt_start = clock();

		//
		//  Step.2 End to Start (EndVertex to LargeDegreeVertexFromS)
		if(Mode_Control==1)  // if Start to End Faild
		{
			//std::cout << "ESmod:::ESmode Step.2 " << std::endl;
			EdgeHashSet* End_set_ES = new EdgeHashSet;
			End_set_ES->set_empty_key("");
			if(LargeDegreeVertexFromS!="")
				End_set_ES->insert(LargeDegreeVertexFromS);    // End_point Hash_set
			else if (LargeDegreeVertexFromS=="")
			{
				std::cout << "mp" << mateNum+1 << "\t Step.2 LargeDegreeVertexFromS is NULL\n" ;
				continue;
			}

			EdgeHashSet* visit_Vertex_set_ES = new EdgeHashSet;
			visit_Vertex_set_ES->set_empty_key("");
			VertexCount=0;
			for(size_t i = 0; i < EndVertixID_list.size(); i++)    // for each start point // left to right search
			{
				if (i>3) break;
				StartVertixID=EndVertixID_list[i];
				pRootVertex = pGraph->getVertex(StartVertixID);	// Get the start(root) vertex from main graph
				curr_distance=0;
				stop=0,	final_distance=0;
				leaf_count=opt::LeafCount-50; 	// no necessary to DFS search , just find LargeDegreeVertexFromE
				//LargeDegreeVertexFromE = "";	// no necessary initial because only keep one time
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
					addNeighborsToSubgraph_dfs(pRootVertex, pSubgraph, StartVertixID_dir, &stop, &final_distance, upper_distance, curr_distance, &leaf_count, LargeDegreeVertexFromE, End_set_ES, visit_Vertex_set_ES, &VertexCount);
				}
				if(stop==1)  // break of muti start
				{
					//std::cout << "mp" << mateNum+1 << "StartVertex: " << StartVertixID << " EndVertex: " << LargeDegreeVertexFromS << " upper_dis= " << upper_distance << std::endl;
					break;
				}
			}  // for start point 

			delete End_set_ES;
			delete visit_Vertex_set_ES;
			//std::cout << "   Step2.VertexCount=" << VertexCount << std::endl;
			VertexCount2=VertexCount;
		}
		if(Mode_Control==1 && final_distance!=0)
			Mode_Control=2; // mode2 = SE failed & ES sucessful
		else if(Mode_Control==1 && final_distance==0)
			Mode_Control=3; // mode3 = SE & ES failed 

		tt_end = clock();
		//std::cout << "=====Time3:Step2:" << ((tt_end-tt_start)/(double)(CLOCKS_PER_SEC)) << std::endl;
		tt_start = clock();
		//
		//  Step.3 Start to End (StartVertex to LargeDegreeVertexFromS)
		std::string outStartVertexID="",outEndVertexID="";
		if(Mode_Control==3)  // if Start to End Faild
		{
			//std::cout << "SEmod:::SEmode Step.3 " << std::endl;
			EdgeHashSet* End_set_SE = new EdgeHashSet;
			End_set_SE->set_empty_key("");
			if(LargeDegreeVertexFromS!="")
				End_set_SE->insert(LargeDegreeVertexFromS);    // End_point Hash_set
			else if (LargeDegreeVertexFromS=="")
			{
				std::cout << "mp" << mateNum+1 << "\t Step.3 LargeDegreeVertexFromS is NULL\n" ;
				continue;
			}

			EdgeHashSet* visit_Vertex_set_SE = new EdgeHashSet;
			visit_Vertex_set_SE->set_empty_key("");
			VertexCount=0;
			//for(size_t i = 0; i < StartVertixID_list.size(); ++i)    // for each start point // Start to End search
			for(size_t i = 0; i < 1; i++)
			{
				//if(i>5) break;
				StartVertixID=StartVertixID_list[i];
				Comp_flagS=0;
				pRootVertex = pGraph->getVertex(StartVertixID);	// Get the start(root) vertex from main graph
				curr_distance=0;
				stop=0,	final_distance=0;
				leaf_count=opt::LeafCount-20;
				upper_distance=opt::depth;

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
					addNeighborsToSubgraph_dfs(pRootVertex, pSubgraph, StartVertixID_dir, &stop, &final_distance, upper_distance, curr_distance, &leaf_count, LargeDegreeVertexFromS, End_set_SE, visit_Vertex_set_SE, &VertexCount);
				}
				if(stop==1)  // break of muti start
				{
					//std::cout << "Step.3 mp" << mateNum+1 << " Flag: "<< Comp_flagS << " StartVertex: " << StartVertixID << " EndVertex: " << LargeDegreeVertexFromS << std::endl;
					outStartVertexID=StartVertixID;
					break;
				}
			}  // for start point 

			delete End_set_SE;
			delete visit_Vertex_set_SE;
			//std::cout << "   Step3.VertexCount=" << VertexCount << std::endl;
			VertexCount3=VertexCount;
		}

		tt_end = clock();
		//std::cout << "=====Time4:Step3:" << ((tt_end-tt_start)/(double)(CLOCKS_PER_SEC)) << std::endl;
		tt_start = clock();

		//
		//  Step.4 End to Start (EndVertex to LargeDegreeVertexFromL)
		//
		if(Mode_Control==3)  // if Start to End Faild
		{
			//std::cout << "ESmod:::ESmode Step.4 " << std::endl;
			EdgeHashSet* End_set_ES = new EdgeHashSet;
			End_set_ES->set_empty_key("");
			if(LargeDegreeVertexFromE!="")
				End_set_ES->insert(LargeDegreeVertexFromE);    // End_point Hash_set
			else if (LargeDegreeVertexFromE=="")
			{
				std::cout << "mp" << mateNum+1 << "\t Step.4 LargeDegreeVertexFromE is NULL\n" ;
				continue;
			}

			EdgeHashSet* visit_Vertex_set_ES = new EdgeHashSet;
			visit_Vertex_set_ES->set_empty_key("");
			VertexCount=0;            
			//for(size_t i = 0; i < EndVertixID_list.size(); ++i)    // for each start point // left to right search
			for(size_t i = 0; i < 1; i++)    // for each start point // left to right search
			{
				//if(i>5) break;
				StartVertixID=EndVertixID_list[i];  // Get the start(root) vertex from main graph
				Comp_flagE=0;
				pRootVertex = pGraph->getVertex(StartVertixID);
				curr_distance=0;
				stop=0,	final_distance=0;
				leaf_count=opt::LeafCount-20;
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
					addNeighborsToSubgraph_dfs(pRootVertex, pSubgraph, StartVertixID_dir, &stop, &final_distance, upper_distance, curr_distance, &leaf_count, LargeDegreeVertexFromE, End_set_ES, visit_Vertex_set_ES, &VertexCount);
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
			//std::cout << "   Step4.VertexCount=" << VertexCount << std::endl;
			VertexCount4=VertexCount;
		}

		tt_end = clock();
		//std::cout << "=====Time5:Step4:" << ((tt_end-tt_start)/(double)(CLOCKS_PER_SEC)) << std::endl;
		tt_start = clock();

		//
		// OutPut: simplify contigs_graph and output contigs
		int thread_num=omp_get_thread_num();
		OutputContig(pSubgraph, Mode_Control, mateNum, final_distance, thread_num, outStartVertexID, outEndVertexID, Comp_flagS, Comp_flagE);

		tt_end = clock();
		//std::cout << "\n=====Time6:OutFile:" << ((tt_end-tt_start)/(double)(CLOCKS_PER_SEC)) << std::endl;
		tt_start = clock();

		delete visit_Vertex_set;
		delete End_set;
		delete pSubgraph;

		tt_end = clock();
		//std::cout << "=====Time7:delete:" << ((tt_end-tt_start)/(double)(CLOCKS_PER_SEC)) << std::endl;

		m_end = clock();
		std::cout << "\tTotal Time:" << ((m_end-m_start)/(double)(CLOCKS_PER_SEC)) ; //<< std::endl; 
		std::cout << " 1:" << VertexCount1 << " 2:" << VertexCount2 << " 3:" << VertexCount3 << " 4:" << VertexCount4 << std::endl;

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
void addNeighborsToSubgraph_dfs(Vertex* pCurrVertex, StringGraph* pSubgraph, EdgeDir dir, int *stop, unsigned int *final_distance, unsigned int upper_distance, unsigned int curr_distance, unsigned int *leaf_count, std::string &LargeDegreeVertex, EdgeHashSet* End_set, EdgeHashSet* visit_Vertex_set, unsigned int *VertexCount)
{

	//clock_t tt_start, tt_end;
	//tt_start = clock();	
	//tt_end = clock();
	//std::cout << "   ===Time:" << ((tt_end-tt_start)/(double)(CLOCKS_PER_SEC)) << std::endl;

	// These are the edges in the main graph
	EdgePtrVec edges=pCurrVertex->getEdges(dir); 
	//std::cout << " vertex:" << pCurrVertex->getID()  << " curr_dis:" << curr_distance << " edge_num:" << edges.size() << std::endl;

	*VertexCount=*VertexCount+1;

	//==================================================================================
	// search control secion
	//==================================================================================

	if ( LargeDegreeVertex == "" && edges.size() > 5 )   // Keep LargeDegreeVertex
	{
		LargeDegreeVertex = pCurrVertex->getID();
		//std::cout << "LargeDegreeVertex: " << LargeDegreeVertex << std::endl;
	}
	/*
	   if(opt::VertexCount > 0 && *VertexCount > opt::VertexCount)  // search control (Vertex number)
	   {
	//std::cout << "TouchVertexCount: " << *VertexCount << std::endl;
	//return;
	//std::cout << "check " << std::endl;
	}	
	*/	

	if(curr_distance > upper_distance)  // search control (depth)
	{
		*leaf_count=*leaf_count+1;
		//std::cout << "leaf_count(depth): " << *leaf_count << std::endl;
		return;
	}	

	if(opt::LeafCount > 0 && *leaf_count > opt::LeafCount) // search control (leaf)
	{
		//std::cout << "out of leaf_count(up_bound): " << *leaf_count << std::endl;
		return;
	}

	if(edges.size()==0)  // if vertex degree=0 leaf_count++
	{
		*leaf_count=*leaf_count+1;
		//std::cout << "leaf_count(degree): " << *leaf_count << std::endl;
		return;
	}

	//tt_end = clock();
	//std::cout << "   ===Time:" << ((tt_end-tt_start)/(double)(CLOCKS_PER_SEC)) << std::endl;
	//tt_start = clock();

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
				*stop=1;
				*final_distance=curr_distance;
				copyVertexToSubgraph(pSubgraph, pCurrVertex);
				copyVertexToSubgraph(pSubgraph, pY);
				Overlap ovr = edges[i]->getOverlap();
				SGAlgorithms::createEdgesFromOverlap(pSubgraph, ovr, true);
				return;
			}

			// Correct dir for complement edge
			assert(edges[i]->getDir() == dir);
			EdgeDir corrected_dir = correctDir(edges[i]->getDir(), edges[i]->getComp()); 
			addNeighborsToSubgraph_dfs(pY, pSubgraph, corrected_dir, stop, final_distance, upper_distance, curr_distance, leaf_count, LargeDegreeVertex, End_set, visit_Vertex_set, VertexCount);
		}
		else if (EHSiter != visit_Vertex_set->end()) // if the vertex has been travaled leaf_count++;
		{
			*leaf_count=*leaf_count+1;	
			//std::cout << "leaf_count(visited): " << *leaf_count << std::endl;
		}


		// Reached EndPoint :: Create subgraph to output Reads
		if(*stop==1)
		{
			Vertex* pY = edges[i]->getEnd();
			copyVertexToSubgraph(pSubgraph, pCurrVertex);
			copyVertexToSubgraph(pSubgraph, pY);
			Overlap ovr = edges[i]->getOverlap();
			SGAlgorithms::createEdgesFromOverlap(pSubgraph, ovr, true);
			return;
		}
	}
}
//
//	output 
//
void OutputContig(StringGraph* pSubgraph, int Mode_Control, unsigned int mateNum, unsigned int final_distance, int thread_num, std::string outStartVertexID, std::string outEndVertexID, int Comp_flagS, int Comp_flagE)
{
	// Simplify and rename SuperReads
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

	// MatePair insert size variance control
	// innerReads ( min_distance <= superRead length <= max_distance)
	// outerReads ( superRead length > max_distance OR superRead length < min_distance)
	unsigned int max_distance = opt::mean_d + opt::var_d ;  // int max_distance = mean + variance ;
	unsigned int min_distance = opt::mean_d - opt::var_d ;  // int min_distance = mean - variance ;

	if(Mode_Control==0 && final_distance<=max_distance && final_distance>=min_distance)
	{	
		// Remove dead-end branches from the graph
		if(opt::numTrimRounds > 0) 
		{
			int numTrims = opt::numTrimRounds;
			while(numTrims-- > 0)
				pSubgraph->visit(trimVisit);
		}

		pSubgraph->simplify_MP(outStartVertexID, outEndVertexID);

		int ContigNum = pSubgraph->simplify_CountVertex();
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
		pSubgraph->simplify_MP(outStartVertexID, outEndVertexID);

		int ContigNum = pSubgraph->simplify_CountVertex();
		if(ContigNum==1)
		{
			ContigName = "mp" + string_mpNum + "_mode0_" ;
			pSubgraph->renameVertices_MP(ContigName, outStartVertexID, outEndVertexID, Comp_flagS, Comp_flagE);
			SGMPFastaVisitor av(outOuterReads);
			pSubgraph->visit(av);
			std::cout << "mp" << mateNum+1 << "_mode0" << "\t Out of dis. real_dis=" << final_distance ;
		}
	}
	else if ( Mode_Control==2 && final_distance!= 0 )
	{

		ContigName = "mp" + string_mpNum + "_mode2_" ; 
		//pSubgraph->simplify();
		//pSubgraph->renameVertices(ContigName);
		pSubgraph->simplify_MP(outStartVertexID, outEndVertexID);
		pSubgraph->renameVertices_MP(ContigName, outStartVertexID, outEndVertexID, Comp_flagS, Comp_flagE);
		SGMPFastaVisitor av(outMode2Reads);
		pSubgraph->visit(av);
		std::cout << "mp" << mateNum+1 << "_mode2" ;

	}
	else if ( Mode_Control==3 && final_distance!= 0 )
	{
		std::string FlagS,FlagE;
		FlagS=int2str(Comp_flagS);
		FlagE=int2str(Comp_flagE);
		pSubgraph->simplify_MP(outStartVertexID, outEndVertexID);

		int ContigNum = pSubgraph->simplify_CountVertex();
		if(ContigNum==2)
		{
			ContigName = "mp" + string_mpNum + "_mode3_" + FlagS + ":" + FlagE + "_" ;
			pSubgraph->renameVertices_MP(ContigName, outStartVertexID, outEndVertexID, Comp_flagS, Comp_flagE);
			SGMPFastaVisitor av(outMode3Reads);
			pSubgraph->visit(av);
			std::cout << "mp" << mateNum+1 << "_mode3" ;
		}
	}
	else
	{
		std::cout << "mp" << mateNum+1 << "\t Can't find a path " ;
	}
}
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
			case 'V': arg >> opt::VertexCount; break;
			case 'L': arg >> opt::trimLengthThreshold; break;
			case 'X': arg >> opt::numTrimRounds; break;
			case '?': die = true; break;
			case 'v': opt::verbose++; break;
			case OPT_TR: opt::bPerformTR = true; break;
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

} 
//
// Load *.fastq to locate Start and End
//
void load_MpSeqInfo()
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
}
// We do not directly add the vertex pointer to the graph, as each graph
// manages its set of vertices and we do not want to double-free the vertices
void copyVertexToSubgraph(StringGraph* pSubgraph, const Vertex* pVertex)
{
	// Make sure the vertex hasn't been added yet
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