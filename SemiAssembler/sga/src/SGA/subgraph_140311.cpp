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
// yuhan
void addNeighborsToSubgraph_dfs(Vertex* pCurrVertex, StringGraph* pSubgraph, EdgeDir dir, unsigned int *stop, unsigned int *final_distance, unsigned int curr_distance, unsigned int leaf_count, EdgeHashSet* End_set, EdgeHashSet* visit_Vertex_set);
void loadMpID();
//void OutputContig(StringGraph* pSubgraph, unsigned int final_distance, unsigned int kkk);
//void graph_recover(Vertex* pCurrVertex, EdgeDir dir, int curr_distance);

// Original
void copyVertexToSubgraph(StringGraph* pSubgraph, const Vertex* pVertex);
//void addNeighborsToSubgraph(Vertex* pCurrVertex, StringGraph* pSubgraph, int span, EdgeDir dir);

//
// Getopt
//
#define SUBPROGRAM "subgraph"
static const char *SUBGRAPH_VERSION_MESSAGE =
SUBPROGRAM " Version " PACKAGE_VERSION "\n"
"Written by Jared Simpson.\n"
"\n"
"Copyright 2010 Wellcome Trust Sanger Institute\n";

static const char *SUBGRAPH_USAGE_MESSAGE =
"Usage: " PACKAGE_NAME " " SUBPROGRAM " [OPTION] ... mp_ID1 mp_ID2 TargetASGQ ASQGFILE\n"
"Extract the subgraph around the sequence with ID from an asqg file with direction.\n"
"\n"
"  -v, --verbose                        display verbose output\n"
"      --help                           display this help and exit\n"
"      -o, --out=FILE                   write the subgraph to FILE (default: subgraph.asqg.gz, subgraph2.asqg.gz)\n"
"      -s, --size=N                     the size of the subgraph to extract, all vertices that are at most N hops\n"
"                                       away from the root will be included (default: 5)\n"
"\nTrimming parameters:\n"
"      -b, --bubble=N                   perform N bubble removal steps (default: 3)\n"
"      -x, --cut-terminal=N             cut off terminal branches in N rounds (default: 10)\n"
"      -l, --min-branch-length=LEN      remove terminal branches only if they are less than LEN bases in length (default: 150)\n"
"      -M, --max-overlap=LEN            keep edges whose overlap length larger than LEN, or keep the longest overlap edge.\n"
"          --transitive-reduction       remove transitive edges from the graph. Off by default.\n"
"\n"
"      -d                               0: sense, 1: antisense, 2: bidirectional\n"
"      -m, --size=N                     mean of insert size (default: 3000)\n"
"      -a, --size=N                     variance of insert size (default: 1000)\n"
"      -t, --size=N                     search depth (default: mean + variance)\n"
"\n"
"      -f, --size=N                     Upper bound of Leaf count (default: 0(all) )\n"
"      -r, --size=N                     Upper bound of Start count (default: 0(all) )\n"
"      -g, --size=N                     Upper bound of MaxEdge of each vertex (default: 0 (all) )\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

namespace opt
{
    static unsigned int verbose;
    static std::string asqgFile;             // Input asqg
    static std::string TasqgFile;            // Input Targer asgq
    static std::string mp_ID1;               // Input matepair ID (xxx/1)
    static std::string mp_ID2;               // Input matepair ID (xxx/2)

    static std::string prefix;               // Output prefix
    static std::string outSubgraph;          // Output Subgraph
    static std::string outContigsFile;       // Output contigs file
    static std::string outGraphFile;         // Output contigs asqg graph
  
    // search path parameters
    static unsigned int Edir = 0;              // 0: sense, 1: antisense, 2: bidirectional
    static unsigned int span = 5;              // subgrph depth
    static unsigned int mean_d = 4500;         // mean of insert size
    static unsigned int var_d = 1500;          // variance of insert size

    // Search Control parameters
    static unsigned int depth = mean_d + var_d + 500;   // subgraph depth
    static unsigned int LeafCount = 0;
    static unsigned int StartCount = 0;       // test
    static unsigned int MaxEdge = 0;          // test

    // Bubble parameters
    //static int numBubbleRounds = 3;

    // Trim parameters
    static int numTrimRounds = 5;
    static unsigned int trimLengthThreshold = 300;
    static unsigned int MaxOverlap = 0;          // test

    // Visitor Control
    static bool bPerformTR  = false;            // SGTransitiveReductionVisitor

    // Data array
    static std::vector<std::string> IDStr1;    // mp_ID1_array
    static std::vector<std::string> IDStr2;    // mp_ID2_array
}

static const char* shortopts = "o:s:d:m:a:x:b:l:i:t:f:r:g:M:";

enum { OPT_HELP = 1, OPT_VERSION, OPT_TR };

static const struct option longopts[] = {
    { "verbose",                  no_argument,       NULL, 'v' },
    { "prefix",                   required_argument, NULL, 'o' },
    //{ "bubble",                   required_argument, NULL, 'b' },
    { "cut-terminal",             required_argument, NULL, 'x' },  // numTrimRounds
    { "min-branch-length",        required_argument, NULL, 'l' },  // trimLengthThreshold
    { "span",                     required_argument, NULL, 's' },  // span
    { "Edir",                     required_argument, NULL, 'd' },  // Edir
    { "mean-of-insert-size",      required_argument, NULL, 'm' },  // mean of insert size
    { "variance-of-insert-size",  required_argument, NULL, 'a' },  // variance of insert size
    { "search-depth",             required_argument, NULL, 't' },  // search depth
    { "leaf-count",               required_argument, NULL, 'f' },  // leaf_count
    { "start-count",              required_argument, NULL, 'r' },  // StartCount
    { "max-edge",                 required_argument, NULL, 'g' },  // MaxEdge
    { "Max-overlap",              required_argument, NULL, 'M' },  // MaxOverlap
    { "transitive-reduction",     no_argument,       NULL, OPT_TR },
    { "help",                     no_argument,       NULL, OPT_HELP },
    { "version",                  no_argument,       NULL, OPT_VERSION },
    { NULL, 0, NULL, 0 }
};

//
// Main
//
int subgraphMain(int argc, char** argv)
{
    Timer* pTimer = new Timer("Total Time");
    parseSubgraphOptions(argc, argv);
    loadMpID(); //load mate pair id
    subgraph();
    delete pTimer;

    return 0; 
}

void subgraph()
{
	// Visitor functors
	SGGraphStatsVisitor statsVisit;
	SGContainRemoveVisitor containVisit;
	SGTransitiveReductionVisitor trVisit;
	SGTrimVisitor trimVisit(opt::trimLengthThreshold);
    SGMaxOverlapVisitor MaxOverlapVisit(opt::MaxOverlap);

    //Timer* pTimer_ADQG = new Timer("loadASQG");
    clock_t t1, t2;
    t1 = clock();
    StringGraph* pGraph = SGUtil::loadASQG(opt::asqgFile, 0, true);
    t2 = clock();
    std::cout << "Done: loadASQG " << "\tTime:" << ((t2-t1)/(double)(CLOCKS_PER_SEC)) << std::endl;
    //delete pTimer_ADQG;

    // Remove containments from the graph
    //std::cout << "Start: SGContainRemoveVisitor" << std::endl;
    //Timer* pTimer_vistor = new Timer("SGContainRemoveVisitor");
    t1 = clock();
	while(pGraph->hasContainment())
	    pGraph->visit(containVisit);
    //delete pTimer_vistor;
    t2 = clock();
    std::cout << "Done: SGContainRemoveVisitor" << "\tTime:" << ((t2-t1)/(double)(CLOCKS_PER_SEC)) << std::endl;
    pGraph->visit(statsVisit); 
    
	// Remove any extraneous transitive edges that may remain in the graph
	if(opt::bPerformTR)
    {
	    t1 = clock();
	    //std::cout << "Start: SGTransitiveReductionVisitor" << std::endl;
	    //pTimer_vistor = new Timer("SGTransitiveReductionVisitor");
	    pGraph->visit(trVisit);
        //delete pTimer_vistor;
        t2 = clock();
        std::cout << "Done: SGTransitiveReductionVisitor" << "\tTime:" << ((t2-t1)/(double)(CLOCKS_PER_SEC)) << std::endl;
        pGraph->visit(statsVisit); 
	}

    //overlap cutoff. if overlap length > maxOverlap, treat it as a repeat.
    if(opt::MaxOverlap > 0)
    {
	    t1 = clock();
	    //std::cout << "Start: SGMaxOverlapVisitor" << std::endl;
	    //pTimer_vistor = new Timer("SGMaxOverlapVisitor");
        pGraph->visit(MaxOverlapVisit);
        //delete pTimer_vistor;
        t2 = clock();
        std::cout << "Done: SGMaxOverlapVisitor" << "\tTime:" << ((t2-t1)/(double)(CLOCKS_PER_SEC)) << std::endl;
        pGraph->visit(statsVisit); 	
    }

    /*
	// Pre-assembly graph stats
	pGraph->visit(statsVisit); 

    //pGraph->printMemSize();
	//pGraph->writeASQG(opt::outSimplifyGraph);
	*/

	//=======================================================================================
	// load target file (matepair map paired reads)
	// load target file (matepair map paired reads)
	//=======================================================================================
    //Timer* pTimer_TASQG = new Timer("load Target ASQG");
    t1 = clock();
    StringGraph* tarGraph = SGUtil::loadASQG(opt::TasqgFile, 0, true);
    t2 = clock();
    std::cout << "Done: load Target ASQG " << "\tTime:" << ((t2-t1)/(double)(CLOCKS_PER_SEC)) << std::endl;
    //delete pTimer_TASQG;

	//=======================================================================================
	// matepair path visitor (for each matepair)
	// matepair path visitor (for each matepair)
	//=======================================================================================
    //int jjj=1; // outfile prefix
    unsigned int kkk=0;
    unsigned int kkk_up=opt::IDStr1.size();
    #pragma omp parallel for num_threads(60) private(kkk) 
	for(kkk=0 ; kkk < kkk_up ; kkk++)  // for each mate pair
	{	
        //Timer* pTimerE = new Timer("Each matepair"); 
        clock_t t1,t2;
        t1 = clock();

        std::string StartVertixID = opt::IDStr1[kkk];  // start ID in mp_ID1 (xxx/1)
	    std::string EndVertixID   = opt::IDStr2[kkk];  // end   ID in mp_ID2 (xxx/2)
        //std::cout << "StartVertixID " << StartVertixID <<  " EndVertixID   " << EndVertixID << std::endl;

        // Get the root vertex (end)
        Vertex* tempVertex_end = tarGraph->getVertex(EndVertixID);
        EdgePtrVec edges_endPoint      = tempVertex_end->getEdges(ED_SENSE);

        EdgeHashSet* End_set = new EdgeHashSet;
        End_set->set_empty_key("");
        for(size_t i = 0; i < edges_endPoint.size(); ++i)
        {
            EndVertixID=(edges_endPoint[i])->getEnd()->getID();
            End_set->insert(EndVertixID.c_str());    // End_point Hash_set
        }
		if(edges_endPoint.size() == 0)
		{
		    std::cout << "mp" << kkk+1 << "\t End Vertex not found.\n";
	        continue;
		}

        // Get the root vertex (start)
	    Vertex* tempRootVertex = tarGraph->getVertex(StartVertixID);
        EdgePtrVec edges_start=tempRootVertex->getEdges(ED_SENSE);  // These are the edges in the target graph
		if(edges_start.size() == 0)
		{
		    std::cout << "mp" << kkk+1 << "\t Start Vertex not found.\n";
	        continue;
		}

		StringGraph* pSubgraph = new StringGraph;
		// Set the graph parameters to match the main graph
		pSubgraph->setContainmentFlag(pGraph->hasContainment());
		pSubgraph->setTransitiveFlag(pGraph->hasTransitive());
		pSubgraph->setMinOverlap(pGraph->getMinOverlap());
		pSubgraph->setErrorRate(pGraph->getErrorRate());

        EdgeHashSet* visit_Vertex_set = new EdgeHashSet;
        visit_Vertex_set->set_empty_key("");

		Vertex* pRootVertex = NULL;
		unsigned int stop=0;
		unsigned int final_distance=0;
        unsigned int curr_distance=0;
        unsigned int leaf_count=0;
        //unsigned int start_count=0;

        //if (opt::StartCount == 0)
        //    start_count = edges_start.size();
//std::cout << "check1" << std::endl;
        for(size_t i = 0; i < edges_start.size(); ++i)    // for each start point
        {
		    if(stop==1)  // break of muti start
                break;
            //if(i > start_count)
            //    break;

		    // Get the root vertex from main graph
	        StartVertixID=(edges_start[i])->getEnd()->getID();
		    pRootVertex = pGraph->getVertex(StartVertixID);

		    stop=0;
		    final_distance=0;
	        curr_distance=0; // initail start distance
		    if(pRootVertex == NULL)
		    {
		        std::cout << "mp" << kkk+1 << "\t Vertex " << StartVertixID << " not found in the graph.\n";
	            continue;
		    }
		    else
		    {
	            EdgeDir StartVertixID_dir = correctDir(edges_start[i]->getDir(), edges_start[i]->getComp());  // correctDir by MP mapping PE

		        // Recursively add neighbors by distance
                addNeighborsToSubgraph_dfs(pRootVertex, pSubgraph, StartVertixID_dir, &stop, &final_distance, curr_distance, leaf_count, End_set, visit_Vertex_set);
	
		        // Write the subgraph
		        //pSubgraph->writeASQG(opt::outSubgraph);
		    }
        }  // for start point 
//std::cout << "check2" << std::endl;
	    // OutPut: simplify contigs_graph and output contigs
        //OutputContig(pSubgraph, final_distance, kkk);

		char bb[5],cc[50];
		sprintf(bb,"%d",omp_get_thread_num());
		std::string tt = bb ; 	    
	    //opt::outContigsFile   = "mp-contig.fa";
	    //opt::outContigsFile   = "mp" + tt + "-contig.fa";
	    std::string outContigsFile = "mp" + tt + "-contig.fa";
        std::string outLargerFile = "mp" + tt + "-outerDistance.fa";

		sprintf(cc,"%d",kkk+1);
		std::string ttt = cc ; 
        std::string ContigName = "mp" + ttt + "_";
	
	    unsigned int max_distance = opt::mean_d + opt::var_d ;  // int max_distance = mean + variance ;
		unsigned int min_distance = opt::mean_d - opt::var_d ;  // int min_distance = mean - variance ;
			
		if(final_distance<=max_distance && final_distance>=min_distance)
		{		
		    // Remove containments from the graph
		    //while(pSubgraph->hasContainment())
		    //    pSubgraph->visit(containVisit);
		
		    // Remove any extraneous transitive edges that may remain in the graph
		    //if(opt::bPerformTR)
		    //    pSubgraph->visit(trVisit);
		    //pSubgraph->visit(trVisit);
		    
		    // Remove dead-end branches from the graph
		    if(opt::numTrimRounds > 0) 
		    {
		        int numTrims = opt::numTrimRounds;
		        while(numTrims-- > 0)
		            pSubgraph->visit(trimVisit);
		        //pSubgraph->visit(statsVisit);
		    }

		    pSubgraph->simplify();
		    pSubgraph->renameVertices(ContigName);   
		    SGMPFastaVisitor av(outContigsFile);
		    pSubgraph->visit(av);
            t2 = clock();
		    std::cout << "mp" << kkk+1 << "\t Find a path. final_distance= " << final_distance << "\tTime:" << ((t2-t1)/(double)(CLOCKS_PER_SEC)) << std::endl;
		}
		else if (final_distance!= 0 && ( final_distance > max_distance || final_distance < min_distance ) )
        {
			pSubgraph->simplify();
		    pSubgraph->renameVertices(ContigName); 
		    SGMPFastaVisitor av(outLargerFile);
		    pSubgraph->visit(av);
            t2 = clock();
		    std::cout << "mp" << kkk+1 << "\t Out of distance. real_distance=" << final_distance << "\tTime:" << ((t2-t1)/(double)(CLOCKS_PER_SEC)) << std::endl;
		}

		else
        {
	        t2 = clock();
		    std::cout << "mp" << kkk+1 << "\t Can't find a path " << "\tTime:" << ((t2-t1)/(double)(CLOCKS_PER_SEC)) << std::endl;
        }

        delete visit_Vertex_set;
        delete End_set;
	    delete pSubgraph;
        //delete pTimerE;
	} // for kkk (each mate pair)
    delete tarGraph;
    delete pGraph;
}
///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////

//////////////// yuhan start ////////////////////////////////////////////////////////////////
//////////////// yuhan start ////////////////////////////////////////////////////////////////
void addNeighborsToSubgraph_dfs(Vertex* pCurrVertex, StringGraph* pSubgraph, EdgeDir dir, unsigned int *stop, unsigned int *final_distance, unsigned int curr_distance, unsigned int leaf_count, EdgeHashSet* End_set, EdgeHashSet* visit_Vertex_set)
{
	
    if(curr_distance > opt::depth)  // search depth control
    {
	    leaf_count++;
	    return;
    }

    // These are the edges in the main graph
    EdgePtrVec edges=pCurrVertex->getEdges(dir); 

    if(opt::MaxEdge != 0 && edges.size() > opt::MaxEdge)
        return;

    if(edges.size()==0)
        leaf_count++;

    if(opt::LeafCount != 0 && leaf_count > opt::LeafCount)
        return;

    for(size_t i = 0; i < edges.size(); ++i)
    {
        Vertex* pY = edges[i]->getEnd();
        // Edge
        //std::string edgeID;
        //std::string edgeID_twin;
        //edgeID = pCurrVertex->getID() + pY->getID();
        //edgeID_twin = pY->getID() + pCurrVertex->getID();
        //EdgeHashSetIter EHSiter = visit_Vertex_set->find(edgeID);
//std::cout << "check1-1" << std::endl;
        // Vertex
        EdgeHashSetIter EHSiter = visit_Vertex_set->find(pY->getID());
//std::cout << "check1-2" << std::endl;
        if( !(EHSiter != visit_Vertex_set->end()) && *stop!=1 )
        { 
            // Edge
            //visit_Vertex_set->insert(edgeID); 
            //visit_Vertex_set->insert(edgeID_twin);

            // Vertex
            visit_Vertex_set->insert(pY->getID());


            // Calculated distance
            curr_distance = curr_distance + (edges[i])->getEnd()->getSeqLen() - (edges[i])->getOverlap().getOverlapLength(0);
            //std::cout << " vertex: " << pY->getID()  << " curr_distance= " << curr_distance << std::endl;
//std::cout << "check1-3" << std::endl;
            // End point hash 
            EdgeHashSetIter RHSiter = End_set->find((edges[i])->getEnd()->getID()); // next vertexi ID
//std::cout << "check1-4" << std::endl;
            if (RHSiter != End_set->end())   // find successful
	        {
                *stop=1;
                *final_distance=curr_distance;
                copyVertexToSubgraph(pSubgraph, pCurrVertex);
                copyVertexToSubgraph(pSubgraph, pY);
                Overlap ovr = edges[i]->getOverlap();
                SGAlgorithms::createEdgesFromOverlap(pSubgraph, ovr, true);
                return;
	        }
//std::cout << "check1-5" << std::endl;
            // search depth control
            //if(curr_distance > opt::depth)   
	        //    return;
            //else
            // Correct dir for complement edge
	    	assert(edges[i]->getDir() == dir);
            EdgeDir corrected_dir = correctDir(edges[i]->getDir(), edges[i]->getComp()); 
            addNeighborsToSubgraph_dfs(pY, pSubgraph, corrected_dir, stop, final_distance, curr_distance, leaf_count, End_set, visit_Vertex_set);
        }
//std::cout << "check1-6" << std::endl;
        // Reached EndPoint
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
/*
void OutputContig(StringGraph* pSubgraph, unsigned int final_distance, unsigned int kkk)
{
	SGGraphStatsVisitor statsVisit;
	SGContainRemoveVisitor containVisit;
	SGTransitiveReductionVisitor trVisit;
	SGTrimVisitor trimVisit(opt::trimLengthThreshold);

	char bb[5];
	sprintf(bb,"%d",omp_get_thread_num());
	std::string tt = bb ; 	    
	//opt::outContigsFile   = "mp-contig.fa";
	//opt::outContigsFile   = "mp" + tt + "-contig.fa";
	std::string outContigsFile = "mp" + tt + "-contig.fa";
    std::string outLargerFile = "mp" + tt + "-outerDistance.fa";
	
	unsigned int max_distance = opt::mean_d + opt::var_d ;  // int max_distance = mean + variance ;
	unsigned int min_distance = opt::mean_d - opt::var_d ;  // int min_distance = mean - variance ;
			
	if(final_distance<=max_distance && final_distance>=min_distance)
	{		
	    // Remove containments from the graph
	    //while(pSubgraph->hasContainment())
	    //    pSubgraph->visit(containVisit);
		
	    // Remove any extraneous transitive edges that may remain in the graph
	    //if(opt::bPerformTR)
	    //    pSubgraph->visit(trVisit);
	    //pSubgraph->visit(trVisit);
		    
	    // Remove dead-end branches from the graph
	    if(opt::numTrimRounds > 0) 
	    {
	        int numTrims = opt::numTrimRounds;
	        while(numTrims-- > 0)
	            pSubgraph->visit(trimVisit);
	        //pSubgraph->visit(statsVisit);
	    }
			    
	    pSubgraph->simplify();
	    pSubgraph->renameVertices("contig-");   
	    SGMPFastaVisitor av(outContigsFile);
	    pSubgraph->visit(av);
	    std::cout << "mp" << kkk+1 << "\t Find a path"<< std::endl;
	}
	else if (final_distance!= 0 && ( final_distance > max_distance || final_distance < min_distance ) )
    {
		pSubgraph->simplify();
	    pSubgraph->renameVertices("OutContig-");   
	    SGMPFastaVisitor av(outLargerFile);
	    pSubgraph->visit(av);
	    std::cout << "mp" << kkk+1 << "\t Out of distance. real_distance=" << final_distance << std::endl; 
	}
	else
	    std::cout << "mp" << kkk+1 << "\t Can't find a path " << std::endl;
}
*/
void loadMpID()
{
	std::string str;
	std::string file1 = opt::mp_ID1;
    std::ifstream fin(file1.c_str(),std::ios::in);
    if (fin == NULL)
    {
	    std::cout << " Can't open mp_ID1 " << std::endl;
        exit(0);
    }
    while (getline(fin, str))
        opt::IDStr1.push_back(str);
    fin.close();

	std::string file2 = opt::mp_ID2;
    std::ifstream fin2(file2.c_str(),std::ios::in);
    if (fin2 == NULL)
    {
	    std::cout << " Can't open mp_ID2 " << std::endl;
        exit(0);
    }
    while (getline(fin2, str))
        opt::IDStr2.push_back(str);
    fin2.close();
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
            case 'l': arg >> opt::trimLengthThreshold; break;
            //case 'b': arg >> opt::numBubbleRounds; break;
            case 'x': arg >> opt::numTrimRounds; break;
            case 'm': arg >> opt::mean_d; break;
            case 'a': arg >> opt::var_d; break;
            case 't': arg >> opt::depth; break;
            case 'f': arg >> opt::LeafCount; break;
            case 'r': arg >> opt::StartCount; break;
            case 'g': arg >> opt::MaxEdge; break;
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
    opt::outGraphFile     = prefix + "-contig-graph.asqg";

    if (argc - optind < 4) 
    {
        std::cerr << SUBPROGRAM ": missing arguments\n";
        die = true;
    } 
    else if (argc - optind > 4) 
    {
        std::cerr << SUBPROGRAM ": too many arguments\n";
        die = true;
    }

    if (die) 
    {
        std::cout << "\n" << SUBGRAPH_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }

    // Parse the vertex id
    opt::mp_ID1 = argv[optind++];
    opt::mp_ID2 = argv[optind++];

    if(opt::mp_ID1.empty())
    {
        std::cerr << SUBPROGRAM ": missing mp_ID1\n";
        exit(EXIT_FAILURE);
    }
    if(opt::mp_ID2.empty())
    {
        std::cerr << SUBPROGRAM ": missing mp_ID2\n";
        exit(EXIT_FAILURE);
    }
    // Parse the input filename
    opt::TasqgFile = argv[optind++];
    opt::asqgFile  = argv[optind++];
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

/*
//
// Original
//
void addNeighborsToSubgraph(Vertex* pCurrVertex, StringGraph* pSubgraph, int span, EdgeDir dir)
{
    if(span <= 0)
        return;

    // These are the edges in the main graph
    EdgePtrVec edges=pCurrVertex->getEdges(dir);
    for(size_t i = 0; i < edges.size(); ++i)
    {    	
        if(edges[i]->getColor() != GC_BLACK)
        {
            Vertex* pY = edges[i]->getEnd();
            
            copyVertexToSubgraph(pSubgraph, pY);
            Overlap ovr = edges[i]->getOverlap();
            SGAlgorithms::createEdgesFromOverlap(pSubgraph, ovr, true);
            edges[i]->setColor(GC_BLACK);
            edges[i]->getTwin()->setColor(GC_BLACK);

            //correct dir for complement edge
	    	assert(edges[i]->getDir() == dir);
            EdgeDir corrected_dir = correctDir(edges[i]->getDir(), edges[i]->getComp());

	        // Recurse
            addNeighborsToSubgraph(pY, pSubgraph, span - 1, corrected_dir);
        }
    }
}
void graph_recover(Vertex* pCurrVertex, EdgeDir dir, int curr_distance)   // make black vertex into white
{
    if(curr_distance > opt::depth)  // search depth control
        return;

    // These are the edges in the main graph
    EdgePtrVec edges=pCurrVertex->getEdges(dir);

    for(size_t i = 0; i < edges.size(); ++i)
    {
        Vertex* pY = edges[i]->getEnd();	

        //if(edges[i]->getColor() == GC_BLACK )
        if(pY->getColor() == GC_BLACK ) 
        {
	        pY->setColor(GC_WHITE);
            //Vertex* pY = edges[i]->getEnd();       
            //edges[i]->setColor(GC_WHITE);
            //edges[i]->getTwin()->setColor(GC_WHITE);
            
            // correct dir for complement edge
	    	assert(edges[i]->getDir() == dir);
            EdgeDir corrected_dir = correctDir(edges[i]->getDir(), edges[i]->getComp());

            // Calculated distance
            curr_distance = curr_distance + (edges[i])->getEnd()->getSeqLen() - (edges[i])->getOverlap().getOverlapLength(0);

            if(curr_distance > opt::depth)  // search depth control
                return;
            else
                graph_recover(pY, corrected_dir, curr_distance);  // Recurse
        }
    }
}
*/
