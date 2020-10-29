#ifndef Settings_H
#define Settings_H

#include <string>


class Settings {

	public:
                        Settings(int argc, char* argv[]);
        double          getAsrvLambda(void) { return asrvLambda; }
        double          getBrlenLambda(void) { return brlenLambda; }
        int             getChainLength(void) { return chainLength; }
        std::string     getDataFilePathName(void) { return dataFilePathName; }
        double          getExpNumCats(void) { return expNumCats; }
        bool            getIsConcFixed(void) { return isConcFixed; }
        std::string     getLogFileName(void) { return logFileName; }
        int             getNumGammaCats(void) { return numGammaCats; }
        std::string     getOutPutFileName(void) { return outPutFileName; }
        int             getPrintFrequency(void) { return printFrequency; }
        double          getPriorConcMean(void) { return priorConcMean; }
        double          getPriorConcVariance(void) { return priorConcVariance; }
        int             getSampleFrequency(void) { return sampleFrequency; }
        std::string     getTreeFileName(void) { return treeFileName; }
        double          getTuningParm(std::string parmNameStr);
        void            setAsrvLambda(double x) { asrvLambda = x; }
        void            setBrlenLambda(double x) { brlenLambda = x; }
        void            setDataFilePathName(std::string s) { dataFilePathName = s; }
        void            setNumGammaCats(int x) { numGammaCats = x; }
        void            setOutPutFileName(std::string s) { outPutFileName = s; }
        void            setTreeFileName(std::string s ) { treeFileName = s; }

    private:
        void            printUsage(void);
        double          brlenLambda;
        int             chainLength;
        int             burninLength;
        std::string     dataFilePathName;
        int             numGammaCats;
        std::string     outPutFileName;
        std::string     treeFileName;
        std::string     logFileName;
        int             printFrequency;
        int             sampleFrequency;
        double          asrvLambda;
        bool            isConcFixed;
        double          expNumCats;
        double          priorConcMean;
        double          priorConcVariance;
        double          tuningParm[5];
};

#endif
