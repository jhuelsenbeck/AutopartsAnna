#include <iostream>
#include <iomanip>
#include <map>
#include <set>
#include <xmmintrin.h>
#include <emmintrin.h>
#include <pmmintrin.h>
#ifdef AP_MPI
#   include <mpi.h>
#endif
#include "Alignment.h"
#include "Chunk.h"
#include "MbMatrix.h"
#include "MbTransitionMatrix.h"
#include "Model.h"
#include "ParmAsrv.h"
#include "ParmFreqs.h"
#include "ParmLength.h"
#include "ParmSubrates.h"
#include "ParmTree.h"
#include "Settings.h"




Chunk::Chunk(Alignment* ap, Settings* sp, Model* mp, int si) {

    pid = 0;
    numProcesses = 1;
#   ifdef AP_MPI
    pid = MPI::COMM_WORLD.Get_rank();
    numProcesses = MPI::COMM_WORLD.Get_size ();
#   endif
    
	// remember the location of important objects
	alignmentPtr = ap;
	modelPtr     = mp;
	settingsPtr  = sp;
	
	// initialize variables
	subsetId     = si;
    update       = false;
	numNodes     = 2 * alignmentPtr->getNumTaxa() - 2;
	numGammaCats = settingsPtr->getNumGammaCats();
	
	// get the list of sites that are included in this "chunk", or subset
	for (int i=0; i<alignmentPtr->getNumChar(); i++)
		{
		if ( alignmentPtr->getIsExcluded(i) == false && alignmentPtr->getPartitionId(i) == subsetId + 1 )
			includedSites.insert(i);
		}
    
    compress();
	

	// allocate the transition probability matrices
	tis = new MbMatrix<double>*[numNodes];
	tis[0] = new MbMatrix<double>[numGammaCats*numNodes];
	for (int i=1; i<numNodes; i++)
		tis[i] = tis[i-1] + numGammaCats;
	for (int i=0; i<numNodes; i++)
		for (int k=0; k<numGammaCats; k++)
			tis[i][k] = MbMatrix<double>(4,4);

	// set up the transition probability calculator
	std::vector<double> bf(4, 0.25);
	std::vector<double> sr(6, 1.0);
	tiMatrix = new MbTransitionMatrix(sr, bf,  true);
		
#	if 0
	std::cout << "Sites for subset " << subsetId << ": ";
	for (std::set<int>::iterator p=includedSites.begin(); p != includedSites.end(); p++)
		std::cout << (*p) << " ";
	std::cout << std::endl;
	printTipCls();
#	endif
}

Chunk::~Chunk(void) {

	delete [] cls;
	delete [] clsPtr;
	delete [] tis[0];
	delete [] tis;
	delete tiMatrix;
}


void Chunk::compress( void ) {

    bool doCompression = true;
    bool treatAmbiguousAsGaps = false;
    
    charMatrix.clear();
    gapMatrix.clear();
    patternCounts.clear();
    numPatterns = 0;
    
    // resize the matrices
    charMatrix.resize(alignmentPtr->getNumTaxa());
    gapMatrix.resize(alignmentPtr->getNumTaxa());
    
    
    // check whether there are ambiguous characters (besides gaps)
    bool ambiguousCharacters = false;

    // find the unique site patterns and compute their respective frequencies
    for (std::set<int>::iterator it = includedSites.begin(); it != includedSites.end(); ++it)
        {
        int site = *it;
        for (int i=0; i<alignmentPtr->getNumTaxa(); i++)
            {
            int nucCode = alignmentPtr->getNucleotide(i, site);
            
            // if we treat unknown characters as gaps and this is an unknown character then we change it
            // because we might then have a pattern more
            if ( treatAmbiguousAsGaps && isAmbiguous(nucCode) )
                {
                nucCode = 15;
                }
            else if ( nucCode != 15 && isAmbiguous(nucCode) )
                {
                ambiguousCharacters = true;
                break;
                }
            }
        
        // break the loop if there was an ambiguous character
        if ( ambiguousCharacters )
            {
            break;
            }
        }
    
    // set the global variable if we use ambiguous characters
    usingAmbiguousCharacters = ambiguousCharacters;
    
    numSites = (int)includedSites.size();
    std::vector<bool> unique(numSites, true);
    // compress the character matrix
    if ( doCompression == true )
        {
        // find the unique site patterns and compute their respective frequencies
        std::map<std::string,size_t> patterns;
        std::set<int>::iterator it = includedSites.begin();
        for (size_t site = 0; site < numSites; ++site, ++it)
            {
            // create the site pattern
            std::string pattern = "";
            for (int i=0; i<alignmentPtr->getNumTaxa(); i++)
                {
                int nucCode = alignmentPtr->getNucleotide(i, *it);
                pattern += nucAsStringValue(nucCode);
                }
            // check if we have already seen this site pattern
            std::map<std::string, size_t>::const_iterator index = patterns.find( pattern );
            if ( index != patterns.end() )
                {
                // we have already seen this pattern
                // increase the frequency counter
                patternCounts[ index->second ]++;
            
                // obviously this site isn't unique nor the first encounter
                unique[site] = false;
                }
            else
                {
                // create a new pattern frequency counter for this pattern
                patternCounts.push_back(1);
                
                // insert this pattern with the corresponding index in the map
                patterns.insert( std::pair<std::string,size_t>(pattern,numPatterns) );
                
                // increase the pattern counter
                numPatterns++;
                
                // flag that this site is unique (or the first occurence of this pattern)
                unique[site] = true;
                }
        
            }
        }
    else
        {
        // we do not compress
        numPatterns = numSites;
        patternCounts = std::vector<size_t>(numSites,1);
        }
    
    
    // allocate and fill the cells of the matrices
    for (int i=0; i<alignmentPtr->getNumTaxa(); i++)
        {
        // resize the column
        charMatrix[i].resize(numPatterns);
        gapMatrix[i].resize(numPatterns);
        size_t patternIndex = 0;
        std::set<int>::iterator it = includedSites.begin();
        for (size_t site = 0; site < numSites; ++site, ++it)
            {
            // only add this site if it is unique
            if ( unique[site] )
                {
                int nucCode = alignmentPtr->getNucleotide(i, *it);
                gapMatrix[i][patternIndex] = (nucCode == 15);
                    
                if ( ambiguousCharacters )
                    {
                    // we use the actual state
                    charMatrix[i][patternIndex] = nucCode;
                    }
                else
                    {
                    // we use the index of the state
                    size_t index = 0;
                    unsigned long state = nucCode;
                    state >>= 1;
                        
                    while ( state != 0 ) // there are still observed states left
                        {
                        // remove this state from the observed states
                        state >>= 1;
                            
                        // increment the index
                        ++index;
                        } // end-while over all observed states for this character
                        
                    charMatrix[i][patternIndex] = index;
                    }
                    
                // increase the pattern index
                patternIndex++;
                }
            }
        }
}

bool Chunk::isAmbiguous(int nc) {
    
    return !(nc == 1 || nc == 2 || nc == 4 || nc == 8);
}

double Chunk::lnLikelihood(bool storeScore) {

	// update the transition probabilities
	updateTransitionProbabilities();
	
	// get some variables that are used over-and-over again
	int stride = 4 * numGammaCats;
	Tree* treePtr = modelPtr->findTree(subsetId);
    size_t pattern_block_start = size_t(floor( (double(pid)   / numProcesses ) * numPatterns) );
    size_t pattern_block_end   = size_t(floor( (double(pid+1) / numProcesses ) * numPatterns) );
		
	/* pass down tree, filling in conditional likelihoods using the sum-product algorithm
	   (Felsenstein pruning algorithm) */
	double* lnScaler = new double[numPatterns];
	for (int i=0; i<numPatterns; i++)
		lnScaler[i] = 0.0;
    
	for (int n=0; n<treePtr->getNumNodes(); n++)
        {
		Node* p = treePtr->getDownPassNode( n );
		if ( p->getIsLeaf() == false )
            {
			if (p->getAnc()->getAnc() == NULL)
                {
				/* three-way split */
				int lftIdx = p->getLft()->getIndex();
				int rhtIdx = p->getRht()->getIndex();
				int ancIdx = p->getAnc()->getIndex();
				int idx    = p->getIndex();
				double *clL = getClsForNode(lftIdx) + pattern_block_start*stride;
				double *clR = getClsForNode(rhtIdx) + pattern_block_start*stride;
				double *clA = getClsForNode(ancIdx) + pattern_block_start*stride;
				double *clP = getClsForNode(idx   ) + pattern_block_start*stride;
				for (size_t c=pattern_block_start; c<pattern_block_end; c++)
                    {
					for (int k=0; k<numGammaCats; k++)
                        {
                        __m128d clL_01 = _mm_load_pd(clL);
                        __m128d clL_23 = _mm_load_pd(clL+2);
                            
                        __m128d clR_01 = _mm_load_pd(clR);
                        __m128d clR_23 = _mm_load_pd(clR+2);
                            
                        __m128d clA_01 = _mm_load_pd(clA);
                        __m128d clA_23 = _mm_load_pd(clA+2);
                        
                        /* A */
                        
                        __m128d a_tis_l_01 = _mm_load_pd(tis[lftIdx][k][0]);
                        __m128d a_tis_l_23 = _mm_load_pd(tis[lftIdx][k][0]+2);
                            
                        __m128d a_tis_r_01 = _mm_load_pd(tis[rhtIdx][k][0]);
                        __m128d a_tis_r_23 = _mm_load_pd(tis[rhtIdx][k][0]+2);
                            
                        __m128d a_tis_a_01 = _mm_load_pd(tis[   idx][k][0]);
                        __m128d a_tis_a_23 = _mm_load_pd(tis[   idx][k][0]+2);
                            
                        __m128d a_l_p01 = _mm_mul_pd(clL_01,a_tis_l_01);
                        __m128d a_l_p23 = _mm_mul_pd(clL_23,a_tis_l_23);
                            
                        __m128d a_r_p01 = _mm_mul_pd(clR_01,a_tis_r_01);
                        __m128d a_r_p23 = _mm_mul_pd(clR_23,a_tis_r_23);
                            
                        __m128d a_a_p01 = _mm_mul_pd(clA_01,a_tis_a_01);
                        __m128d a_a_p23 = _mm_mul_pd(clA_23,a_tis_a_23);
                            
                        __m128d a_sum_L = _mm_hadd_pd(a_l_p01,a_l_p23);
                        __m128d a_sum_R = _mm_hadd_pd(a_r_p01,a_r_p23);
                        __m128d a_sum_A = _mm_hadd_pd(a_a_p01,a_a_p23);
                            
                        /* C */
                            
                        __m128d c_tis_l_01 = _mm_load_pd(tis[lftIdx][k][1]);
                        __m128d c_tis_l_23 = _mm_load_pd(tis[lftIdx][k][1]+2);
                            
                        __m128d c_tis_r_01 = _mm_load_pd(tis[rhtIdx][k][1]);
                        __m128d c_tis_r_23 = _mm_load_pd(tis[rhtIdx][k][1]+2);
                            
                        __m128d c_tis_a_01 = _mm_load_pd(tis[   idx][k][1]);
                        __m128d c_tis_a_23 = _mm_load_pd(tis[   idx][k][1]+2);
                            
                        __m128d c_l_p01 = _mm_mul_pd(clL_01,c_tis_l_01);
                        __m128d c_l_p23 = _mm_mul_pd(clL_23,c_tis_l_23);
                            
                        __m128d c_r_p01 = _mm_mul_pd(clR_01,c_tis_r_01);
                        __m128d c_r_p23 = _mm_mul_pd(clR_23,c_tis_r_23);
                            
                        __m128d c_a_p01 = _mm_mul_pd(clA_01,c_tis_a_01);
                        __m128d c_a_p23 = _mm_mul_pd(clA_23,c_tis_a_23);
                            
                        __m128d c_sum_L = _mm_hadd_pd(c_l_p01,c_l_p23);
                        __m128d c_sum_R = _mm_hadd_pd(c_r_p01,c_r_p23);
                        __m128d c_sum_A = _mm_hadd_pd(c_a_p01,c_a_p23);
                            
                        /* G */
                            
                        __m128d g_tis_l_01 = _mm_load_pd(tis[lftIdx][k][2]);
                        __m128d g_tis_l_23 = _mm_load_pd(tis[lftIdx][k][2]+2);
                        
                        __m128d g_tis_r_01 = _mm_load_pd(tis[rhtIdx][k][2]);
                        __m128d g_tis_r_23 = _mm_load_pd(tis[rhtIdx][k][2]+2);
                            
                        __m128d g_tis_a_01 = _mm_load_pd(tis[   idx][k][2]);
                        __m128d g_tis_a_23 = _mm_load_pd(tis[   idx][k][2]+2);
                            
                        __m128d g_l_p01 = _mm_mul_pd(clL_01,g_tis_l_01);
                        __m128d g_l_p23 = _mm_mul_pd(clL_23,g_tis_l_23);
                            
                        __m128d g_r_p01 = _mm_mul_pd(clR_01,g_tis_r_01);
                        __m128d g_r_p23 = _mm_mul_pd(clR_23,g_tis_r_23);
                            
                        __m128d g_a_p01 = _mm_mul_pd(clA_01,g_tis_a_01);
                        __m128d g_a_p23 = _mm_mul_pd(clA_23,g_tis_a_23);
                            
                        __m128d g_sum_L = _mm_hadd_pd(g_l_p01,g_l_p23);
                        __m128d g_sum_R = _mm_hadd_pd(g_r_p01,g_r_p23);
                        __m128d g_sum_A = _mm_hadd_pd(g_a_p01,g_a_p23);
                            
                        /* T */
                            
                        __m128d t_tis_l_01 = _mm_load_pd(tis[lftIdx][k][3]);
                        __m128d t_tis_l_23 = _mm_load_pd(tis[lftIdx][k][3]+2);
                            
                        __m128d t_tis_r_01 = _mm_load_pd(tis[rhtIdx][k][3]);
                        __m128d t_tis_r_23 = _mm_load_pd(tis[rhtIdx][k][3]+2);
                            
                        __m128d t_tis_a_01 = _mm_load_pd(tis[   idx][k][3]);
                        __m128d t_tis_a_23 = _mm_load_pd(tis[   idx][k][3]+2);
                            
                        __m128d t_l_p01 = _mm_mul_pd(clL_01,t_tis_l_01);
                        __m128d t_l_p23 = _mm_mul_pd(clL_23,t_tis_l_23);
                            
                        __m128d t_r_p01 = _mm_mul_pd(clR_01,t_tis_r_01);
                        __m128d t_r_p23 = _mm_mul_pd(clR_23,t_tis_r_23);
                            
                        __m128d t_a_p01 = _mm_mul_pd(clA_01,t_tis_a_01);
                        __m128d t_a_p23 = _mm_mul_pd(clA_23,t_tis_a_23);
                            
                        __m128d t_sum_L = _mm_hadd_pd(t_l_p01,t_l_p23);
                        __m128d t_sum_R = _mm_hadd_pd(t_r_p01,t_r_p23);
                        __m128d t_sum_A = _mm_hadd_pd(t_a_p01,t_a_p23);
                            
                            
                        // now put it all together

                        __m128d ac_sum_L = _mm_hadd_pd(a_sum_L,c_sum_L);
                        __m128d ac_sum_R = _mm_hadd_pd(a_sum_R,c_sum_R);
                        __m128d ac_sum_A = _mm_hadd_pd(a_sum_A,c_sum_A);
                        __m128d gt_sum_L = _mm_hadd_pd(g_sum_L,t_sum_L);
                        __m128d gt_sum_R = _mm_hadd_pd(g_sum_R,t_sum_R);
                        __m128d gt_sum_A = _mm_hadd_pd(g_sum_A,t_sum_A);
                            
                            
                        __m128d ac_prod_LR = _mm_mul_pd(ac_sum_L,ac_sum_R);
                        __m128d gt_prod_LR = _mm_mul_pd(gt_sum_L,gt_sum_R);
                            
                        __m128d ac = _mm_mul_pd(ac_prod_LR,ac_sum_A);
                        __m128d gt = _mm_mul_pd(gt_prod_LR,gt_sum_A);

                        _mm_store_pd(clP,ac);
                        _mm_store_pd(clP+2,gt);
                            
                        clP += 4;
                        clL += 4;
                        clR += 4;
                        clA += 4;
                        }
                    }
                }
			else
                {
				/* two-way split */
				int lftIdx = p->getLft()->getIndex();
				int rhtIdx = p->getRht()->getIndex();
				int idx    = p->getIndex();
				double *clL = getClsForNode(lftIdx) + pattern_block_start*stride;
				double *clR = getClsForNode(rhtIdx) + pattern_block_start*stride;
                double *clP = getClsForNode(idx   ) + pattern_block_start*stride;
                for (size_t c=pattern_block_start; c<pattern_block_end; c++)
                    {
					for (int k=0; k<numGammaCats; k++)
                        {
                            
                        __m128d clL_01 = _mm_load_pd(clL);
                        __m128d clL_23 = _mm_load_pd(clL+2);
                        
                        __m128d clR_01 = _mm_load_pd(clR);
                        __m128d clR_23 = _mm_load_pd(clR+2);
                        
                        /* A */
                            
                        __m128d a_tis_l_01 = _mm_load_pd(tis[lftIdx][k][0]);
                        __m128d a_tis_l_23 = _mm_load_pd(tis[lftIdx][k][0]+2);
                            
                        __m128d a_tis_r_01 = _mm_load_pd(tis[rhtIdx][k][0]);
                        __m128d a_tis_r_23 = _mm_load_pd(tis[rhtIdx][k][0]+2);
                            
                        __m128d a_l_p01 = _mm_mul_pd(clL_01,a_tis_l_01);
                        __m128d a_l_p23 = _mm_mul_pd(clL_23,a_tis_l_23);
                            
                        __m128d a_r_p01 = _mm_mul_pd(clR_01,a_tis_r_01);
                        __m128d a_r_p23 = _mm_mul_pd(clR_23,a_tis_r_23);
                            
                        __m128d a_sum_L = _mm_hadd_pd(a_l_p01,a_l_p23);
                        __m128d a_sum_R = _mm_hadd_pd(a_r_p01,a_r_p23);
                            
                        /* C */
                            
                        __m128d c_tis_l_01 = _mm_load_pd(tis[lftIdx][k][1]);
                        __m128d c_tis_l_23 = _mm_load_pd(tis[lftIdx][k][1]+2);
                            
                        __m128d c_tis_r_01 = _mm_load_pd(tis[rhtIdx][k][1]);
                        __m128d c_tis_r_23 = _mm_load_pd(tis[rhtIdx][k][1]+2);
                            
                        __m128d c_l_p01 = _mm_mul_pd(clL_01,c_tis_l_01);
                        __m128d c_l_p23 = _mm_mul_pd(clL_23,c_tis_l_23);
                            
                        __m128d c_r_p01 = _mm_mul_pd(clR_01,c_tis_r_01);
                        __m128d c_r_p23 = _mm_mul_pd(clR_23,c_tis_r_23);
                            
                        __m128d c_sum_L = _mm_hadd_pd(c_l_p01,c_l_p23);
                        __m128d c_sum_R = _mm_hadd_pd(c_r_p01,c_r_p23);
                            
                        /* G */
                            
                        __m128d g_tis_l_01 = _mm_load_pd(tis[lftIdx][k][2]);
                        __m128d g_tis_l_23 = _mm_load_pd(tis[lftIdx][k][2]+2);
                            
                        __m128d g_tis_r_01 = _mm_load_pd(tis[rhtIdx][k][2]);
                        __m128d g_tis_r_23 = _mm_load_pd(tis[rhtIdx][k][2]+2);
                            
                        __m128d g_l_p01 = _mm_mul_pd(clL_01,g_tis_l_01);
                        __m128d g_l_p23 = _mm_mul_pd(clL_23,g_tis_l_23);
                            
                        __m128d g_r_p01 = _mm_mul_pd(clR_01,g_tis_r_01);
                        __m128d g_r_p23 = _mm_mul_pd(clR_23,g_tis_r_23);
                            
                        __m128d g_sum_L = _mm_hadd_pd(g_l_p01,g_l_p23);
                        __m128d g_sum_R = _mm_hadd_pd(g_r_p01,g_r_p23);
                            
                        /* T */
                            
                        __m128d t_tis_l_01 = _mm_load_pd(tis[lftIdx][k][3]);
                        __m128d t_tis_l_23 = _mm_load_pd(tis[lftIdx][k][3]+2);
                            
                        __m128d t_tis_r_01 = _mm_load_pd(tis[rhtIdx][k][3]);
                        __m128d t_tis_r_23 = _mm_load_pd(tis[rhtIdx][k][3]+2);
                            
                        __m128d t_l_p01 = _mm_mul_pd(clL_01,t_tis_l_01);
                        __m128d t_l_p23 = _mm_mul_pd(clL_23,t_tis_l_23);
                            
                        __m128d t_r_p01 = _mm_mul_pd(clR_01,t_tis_r_01);
                        __m128d t_r_p23 = _mm_mul_pd(clR_23,t_tis_r_23);
                            
                        __m128d t_sum_L = _mm_hadd_pd(t_l_p01,t_l_p23);
                        __m128d t_sum_R = _mm_hadd_pd(t_r_p01,t_r_p23);
                            
                            
                        // now put it all together
                            
                        __m128d ac_sum_L = _mm_hadd_pd(a_sum_L,c_sum_L);
                        __m128d ac_sum_R = _mm_hadd_pd(a_sum_R,c_sum_R);
                        __m128d gt_sum_L = _mm_hadd_pd(g_sum_L,t_sum_L);
                        __m128d gt_sum_R = _mm_hadd_pd(g_sum_R,t_sum_R);
                            
                            
                        __m128d ac = _mm_mul_pd(ac_sum_L,ac_sum_R);
                        __m128d gt = _mm_mul_pd(gt_sum_L,gt_sum_R);
                                                        
                        _mm_store_pd(clP,ac);
                        _mm_store_pd(clP+2,gt);
                            
                        clP += 4;
                        clL += 4;
                        clR += 4;
                        
                        }
                    }
                }
            
				
			/* scale */
#			if 1
            double *clP = getClsForNode( p->getIndex() ) + pattern_block_start*stride;
            for (size_t c=pattern_block_start; c<pattern_block_end; c++)
                {
				double maxVal = 0.0;
				for (int i=0; i<stride; i++)
                    {
					if (clP[i] > maxVal)
                        {
						maxVal = clP[i];
                        }
                    }
				double scaler = 1.0 / maxVal;
				for (int i=0; i<stride; i++)
                    {
					clP[i] *= scaler;
                    }
				lnScaler[c] += log(maxVal);
				clP += stride;
                }
#			endif
						
        }
    }
		
	/* calculate likelihood */
	Node* p = treePtr->getRoot()->getLft();
	std::vector<double> f = modelPtr->findBaseFreqs(subsetId)->getFreq();
	double catProb = 1.0 / numGammaCats;
	double *clP = getClsForNode( p->getIndex() ) + pattern_block_start*stride;
    double lnL = 0.0;
    for (size_t c=pattern_block_start; c<pattern_block_end; c++)
        {
		double like = 0.0;
		for (int k=0; k<numGammaCats; k++)
            {
			for (int i=0; i<4; i++)
                {
				like += clP[i] * f[i] * catProb;
                }
            clP += 4;
            }
        lnL += (log( like ) + lnScaler[c]) * patternCounts[c];
        }
		
    delete [] lnScaler;
    
#   ifdef AP_MPI
    // std::cerr << "Before: lnL[pid=" << pid << "] = " << lnL << std::endl;
    
    MPI::COMM_WORLD.Barrier();
    
    if ( pid != 0 )
        {
        // send from the workers the log-likelihood to the master
        MPI::COMM_WORLD.Send(&lnL, 1, MPI::DOUBLE, 0, 0);
        }
    
    if ( pid == 0 )
        {
        for (size_t i=1; i<numProcesses; ++i)
            {
            double tmp = 0;
            MPI::COMM_WORLD.Recv(&tmp, 1, MPI::DOUBLE, (int)i, 0);
            lnL += tmp;
            }
        }
    
    MPI::COMM_WORLD.Barrier();
    // MPI::COMM_WORLD.Bcast(&lnL, 1, MPI::DOUBLE, 0);
    
    if ( pid == 0 )
        {
        for (size_t i=1; i<numProcesses; ++i)
            {
            MPI::COMM_WORLD.Send(&lnL, 1, MPI::DOUBLE, (int)i, 0);
            }
        }
    else
        {
        MPI::COMM_WORLD.Recv(&lnL, 1, MPI::DOUBLE, 0, 0);
        }
    
    MPI::COMM_WORLD.Barrier();
    // std::cerr << "After: lnL[pid=" << pid << "] = " << lnL << std::endl;

#endif
	
    // std::cerr << "lnL[" << subsetId << "] = " << lnL << std::endl;
    // std::cerr << "\nExpected:\n" << "lnL[" << subsetId << "] = " << -1490.98 << std::endl;
    // std::exit(0);

    
    if (storeScore == true)
        storedLnL = lnL;
    update = false;
    mostRecentLnL = lnL;
    
	return lnL;
}


/*-------------------------------------------------------------------
 |
 |   GetPossibleNucs:
 |
 |   This function initializes a vector, nuc[MAX_NUM_STATES]. The four elements
 |   of nuc correspond to the four nucleotides in alphabetical order.
 |   We are assuming that the nucCode is a binary representation of
 |   the nucleotides that are consistent with the observation. For
 |   example, if we observe an A, then the nucCode is 1 and the
 |   function initalizes nuc[0] = 1 and the other elements of nuc
 |   to be 0.
 |
 |   Observation    nucCode        nuc
 |        A            1           1000
 |        C            2           0100
 |        G            4           0010
 |        T            8           0001
 |        R            5           1010
 |        Y           10           0101
 |        M            3           1100
 |        K           12           0011
 |        S            6           0110
 |        W            9           1001
 |        H           11           1101
 |        B           14           0111
 |        V            7           1110
 |        D           13           1011
 |        N - ?       15           1111
 |
 -------------------------------------------------------------------*/
std::string Chunk::nucAsStringValue(int nucCode) {

    std::string nucString = "-";
    if (nucCode == 1)
        {
        nucString = "A";
        }
    else if (nucCode == 2)
        {
        nucString = "C";
        }
    else if (nucCode == 3)
        {
        nucString = "M";
        }
    else if (nucCode == 4)
        {
        nucString = "G";
        }
    else if (nucCode == 5)
        {
        nucString = "R";
        }
    else if (nucCode == 6)
        {
        nucString = "S";
        }
    else if (nucCode == 7)
        {
        nucString = "V";
        }
    else if (nucCode == 8)
        {
        nucString = "T";
        }
    else if (nucCode == 9)
        {
        nucString = "W";
        }
    else if (nucCode == 10)
        {
        nucString = "Y";
        }
    else if (nucCode == 11)
        {
        nucString = "H";
        }
    else if (nucCode == 12)
        {
        nucString = "K";
        }
    else if (nucCode == 13)
        {
        nucString = "D";
        }
    else if (nucCode == 14)
        {
        nucString = "B";
        }
    else if (nucCode == 15)
        {
        nucString = "-";
        }
    return nucString;
}

void Chunk::printTis(void) {

	for (int nde=0; nde<numNodes; nde++)
        {
        std::cout << "Transition probabilities for node " << nde << std::endl;
        for (int i=0; i<4; i++)
            {
            std::cout << "   ";
            for (int k=0; k<numGammaCats; k++)
                {
                for (int j=0; j<4; j++)
                    std::cout << std::fixed << std::setprecision(4) << tis[nde][k][i][j] << " ";
                if (k + 1 != numGammaCats)
                    std::cout << " -- ";
                }
            std::cout << std::endl;
            }
        }
}

void Chunk::printTipCls(void) {

	for (int c=0; c<numSites; c++)
		{
		for (int k=0; k<numGammaCats; k++)
			{
			if (k == 0)
				std::cout << std::setw(4) << c << " -- ";
			else 
				std::cout << "        ";
			for (int i=0; i<alignmentPtr->getNumTaxa(); i++)
				{
				double *clp = clsPtr[i] + (c * 4 * numGammaCats) + 4*k;
				for (int s=0; s<4; s++)
					std::cout << std::fixed << std::setprecision(0) << clp[s];
				std::cout << " ";
				}
			std::cout << std::endl;
			}
		}
}

void Chunk::updateTransitionProbabilities(void) {

	// get poiters to the parameters that are necessary to properly update the transition probabilities
	Tree*       treePtr       = modelPtr->findTree(subsetId);
	TreeLength* treeLengthPtr = modelPtr->findTreeLength(subsetId);
	SubRates*   subRatesPtr   = modelPtr->findSubRates(subsetId);
	BaseFreqs*  baseFreqsPtr  = modelPtr->findBaseFreqs(subsetId);
	Asrv*       asrvPtr       = modelPtr->findAsrv(subsetId);
	
	// get the parameter values
	std::vector<double> bf = baseFreqsPtr->getFreq();
	std::vector<double> sr = subRatesPtr->getSubRate();
	std::vector<double> r  = asrvPtr->getRate();
	double len             = treeLengthPtr->getLength();
    
	// update the rate matrix (also updates the eigensystem
	tiMatrix->updateQ(sr, bf);

	// update the transition probabilities
	for (int n=0; n<treePtr->getNumNodes(); n++)
		{
		Node* p = treePtr->getDownPassNode( n );
		if ( p->getAnc() != NULL )
			{
			int idx = p->getIndex();
			double prop = p->getP();
			double v    = len * prop;
			for (int k=0; k<numGammaCats; k++)
				{
				double vr = v * r[k];
				tis[idx][k] = tiMatrix->tiProbs(vr, tis[idx][k]);
				//std::cout << k << " " << vr << std::endl;
				//cout << ti[idx][k] << endl;
				}
			}
		}
        
#   if 0
    std::cout << "bf: ";
    for (int i=0; i<bf.size(); i++)
        std::cout << bf[i] << " ";
    std::cout << std::endl;
    std::cout << "sr: ";
    for (int i=0; i<sr.size(); i++)
        std::cout << sr[i] << " ";
    std::cout << std::endl;
    std::cout << "r: ";
    for (int i=0; i<r.size(); i++)
        std::cout << r[i] << " ";
    std::cout << std::endl;
    std::cout << "len: " << len << std::endl;
#   endif
}



