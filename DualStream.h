#ifndef DualStream_H
#define DualStream_H

#include <iostream>
#include <ostream>


#ifdef AP_MPI
#   include <mpi.h>
#endif

class DualStream {

    public:
        DualStream(std::ostream& os1, std::ostream& os2) : os1(os1), os2(os2) { }
        template<class T>
        DualStream& operator<<(const T& x)
            {
            size_t pid = 0;
#           ifdef AP_MPI
            pid = MPI::COMM_WORLD.Get_rank();
#           endif
            if ( pid == 0 )
                {
                    os1 << x;
                    os2 << x;
                }
            return *this;
            }
            
    private:
        std::ostream&       os1;
        std::ostream&       os2;
};

#endif
