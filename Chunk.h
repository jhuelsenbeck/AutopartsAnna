#ifndef Chunk_H
#define Chunk_H

#include <set>
#include <vector>
#include "MbMatrix.h"

class Alignment;
class MbTransitionMatrix;
class Model;
class Settings;

class Chunk {

	public:
                                            Chunk(Alignment* ap, Settings* sp, Model* mp, int si);
                                           ~Chunk(void);
        double*                             getClsForNode(int i) { return clsPtr[i]; }
        double                              getStoredLnL(void) { return storedLnL; }
        bool                                getUpdate(void) { return update; }
        double                              lnLikelihood(bool storeScore);
        void                                setStoredLnLToMostRecentLnL(void) { storedLnL = mostRecentLnL; }
        void                                setStoredLnL(double x) { storedLnL = x; }
        void                                setUpdate(bool tf) { update = tf; }
        void                                updateTransitionProbabilities(void);

    private:
        void                                compress(void);
        bool                                isAmbiguous(int nc);
        std::string                         nucAsStringValue(int nc);
        void                                printTipCls(void);
        void                                printTis(void);
        size_t                              pid;
        size_t                              numProcesses;
        Settings*                           settingsPtr;
        Alignment*                          alignmentPtr;
        Model*                              modelPtr;
        std::set<int>                       includedSites;
        int                                 subsetId;
        int                                 numSites;
        int                                 numNodes;
        int                                 numGammaCats;
        double*                             cls;
        double**                            clsPtr;
        MbMatrix<double>**                  tis;
        MbTransitionMatrix*                 tiMatrix;
        double                              storedLnL;
        double                              mostRecentLnL;
        bool                                update;
        bool                                usingAmbiguousCharacters;
        std::vector<std::vector<size_t> >   charMatrix;
        std::vector<std::vector<bool> >     gapMatrix;
        int                                 numPatterns;
        std::vector<size_t>                 patternCounts;
};

#endif
