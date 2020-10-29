#ifndef Restaurant_H
#define Restaurant_H

#include <fstream>
#include <map>
#include <set>
#include "Table.h"

class Alignment;
class DualStream;
class MbRandom;
class Model;
class Parm;
class Settings;

class Restaurant {

	public:
                                Restaurant(MbRandom* rp, Settings* sp, Alignment* ap, Model* mp, DualStream* lg, ParmId pid, int np, double ps);
                               ~Restaurant(void);
        bool                    change(void);
        bool                    changeParmOnTable(Table* tbl);
        std::string             getName(void) { return name; }
        int                     getNumTables(void) { return (int)tables.size(); }
        ParmId                  getParmId(void) { return parmId; }
        Table*                  getTableWithPatron(int patron);
        void                    print(void);
        void                    saveState(int n);
        void                    updateSeating(int patron);

    private:
        double                  acceptanceProb(double lnR);
        double                  calcAlphaFromEt(double expT);
        void                    deleteUnoccupiedTables(void);
        double                  expNumTables(double a);
        void                    normalizeLogProbabilitiesInMap(std::map<size_t,double>& m);
        Table*                  pickTableAtRandomFromPrior(void);
        Table*                  pickTableUniformlyAtRandom(void);
        void                    removeTable(Table* tbl);
        double                  sampleAlpha(int k, int n, double oldAlpha, double a, double b);
        void                    setRgf(void);
        std::vector<Table*>     tables;
        Alignment*              alignmentPtr;
        MbRandom*               ranPtr;
        Model*                  modelPtr;
        Settings*               settingsPtr;
        DualStream*             outLog;
        ParmId                  parmId;
        std::string             name;
        int                     numPatrons;
        double                  alpha;
        double                  gammaAlpha;
        double                  gammaBeta;
        double                  probUpdateSeating;
        int                     numAuxTables;
        std::ofstream           parmOut;
        std::ofstream           partOut;
        int*                    rgf;
};

#endif
