#ifndef Parm_H
#define Parm_H

#include <string>

class DualStream;
class MbRandom;
class Model;

class Parm {

	public:
                                Parm(MbRandom* rp, Model* mp, DualStream* lg, std::string nm);
        virtual                ~Parm(void);
        Parm&                   operator=(Parm& b);
        virtual void            print(void)=0;
        virtual double          update(void)=0;
        virtual double          lnPriorProb(void)=0;
        virtual std::string     getParmString(int n)=0;
        std::string             getParmName(void) { return parmName; }
        virtual std::string     getParmHeader(int n)=0;

    protected:
        MbRandom*               ranPtr;
        Model*                  modelPtr;
        DualStream*             outLog;
        std::string             parmName;
};

#endif
