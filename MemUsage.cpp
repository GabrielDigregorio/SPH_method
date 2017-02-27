#include "SPH.hpp"

#if defined(_WIN32) || defined(WIN32) || defined(__WIN32) && !defined(__CYGWIN__)
#include <windows.h>
#include <psapi.h>
    #include <winbase.h>
    #include <stdio.h>
    #pragma comment( lib, "psapi.lib" )

#elif defined(__unix__) || defined(__unix) || defined(unix) || (defined(__APPLE__) && defined(__MACH__))
#include <unistd.h>
#include <sys/resource.h>
    #include "sys/types.h"
    #include "sys/sysinfo.h"

#if defined(__APPLE__) && defined(__MACH__)
#include <mach/mach.h>

#elif defined(__linux__) || defined(__linux) || defined(linux) || defined(__gnu_linux__)
#include <stdio.h>
    #include "sys/types.h"
    #include "sys/sysinfo.h"
    struct sysinfo memInfo;
    
#endif

#else
#error "Cannot define GetMemory( ) or GetMemoryProcessPeak( ) or GetMemoryProcess() for an unknown OS."
#endif


size_t GetMemory(bool screen, bool print)
{
    // open a file to write the memory consumption
    std::ofstream myfile;
    myfile.open ("Memory.txt", std::ofstream::out | std::ofstream::app);

    #if defined(_WIN32) || defined(WIN32) || defined(__WIN32) && !defined(__CYGWIN__)
        /* Windows -------------------------------------------------- */
        // Memory request
        MEMORYSTATUSEX memInfo;
        memInfo.dwLength = sizeof(MEMORYSTATUSEX);
        GlobalMemoryStatusEx(&memInfo);

            // Total Virtual Memory
            DWORDLONG totalVirtualMem = memInfo.ullTotalPageFile;
            // Virtual Memory currently used
            DWORDLONG virtualMemUsed = memInfo.ullTotalPageFile - memInfo.ullAvailPageFile;
            //Total Physical Memory (RAM)
            DWORDLONG totalPhysMem = memInfo.ullTotalPhys;
            //Physical Memory currently used:
            DWORDLONG physMemUsed = memInfo.ullTotalPhys - memInfo.ullAvailPhys;


        if(screen == true){
            std::cout<<"\n TotVirtMem:    \t"<< totalVirtualMem <<" [B]\n";
            std::cout<<" VirtMemCurrUsed:\t"<< virtualMemUsed <<" [B]\n";
            std::cout<<" TotPhysMem(RAM):\t"<< totalPhysMem <<" [B]\n";
            std::cout<<" PhysMemCurrUsed:\t"<< physMemUsed <<" [B]\n \n";
        }
        if(print == true)    
            myfile << totalVirtualMem << " " << virtualMemUsed << " " << totalPhysMem << " " << physMemUsed << "\n" ;

        myfile.close();
        return (size_t) physMemUsed;

    #elif defined(__unix__) || defined(__unix) || defined(unix)
        /* BSD, Linux, and OSX -------------------------------------- */
        // Memory request
        sysinfo (&memInfo);

            // Total Virtual Memory
            long long totalVirtualMem = memInfo.totalram;
            //Add other values in next statement to avoid int overflow on right hand side...
            totalVirtualMem += memInfo.totalswap;
            totalVirtualMem *= memInfo.mem_unit;

            // Virtual Memory currently used
            long long virtualMemUsed = memInfo.totalram - memInfo.freeram;
            //Add other values in next statement to avoid int overflow on right hand side...
            virtualMemUsed += memInfo.totalswap - memInfo.freeswap;
            virtualMemUsed *= memInfo.mem_unit;

            //Total Physical Memory (RAM)
            long long totalPhysMem = memInfo.totalram;
            //Multiply in next statement to avoid int overflow on right hand side...
            totalPhysMem *= memInfo.mem_unit;

            //Physical Memory currently used:
            long long physMemUsed = memInfo.totalram - memInfo.freeram;
            //Multiply in next statement to avoid int overflow on right hand side...
            physMemUsed *= memInfo.mem_unit;

        if(screen == true){
            std::cout<<"\n TotVirtMem:    \t"<< totalVirtualMem <<" [B]\n";
            std::cout<<" VirtMemCurrUsed:\t"<< virtualMemUsed <<" [B]\n";
            std::cout<<" TotPhysMem(RAM):\t"<< totalPhysMem <<" [B]\n";
            std::cout<<" PhysMemCurrUsed:\t"<< physMemUsed <<" [B]\n \n";
        }
        if(print == true)    
            myfile << totalVirtualMem << " " << virtualMemUsed << " " << totalPhysMem << " " << physMemUsed << "\n" ;

        myfile.close();
        return (size_t) physMemUsed;

    #else
        /* MAC, AIX, BSD, Solaris, and Unknown OS ------------------------ */
        return (size_t)0L;			/* Unsupported. */

    #endif

}


/*
 * Returns the peak (maximum so far) resident set size (physical
 * memory use) measured in bytes, or zero if the value cannot be
 * determined on this OS.
 */
size_t GetMemoryProcessPeak(bool screen, bool print)
{
    // open a file to write the memory consumption
    std::ofstream myfile;
    myfile.open ("MemoryProcessPeak.txt", std::ofstream::out | std::ofstream::app);

    #if defined(_WIN32) || defined(WIN32) || defined(__WIN32) && !defined(__CYGWIN__)
        /* Windows -------------------------------------------------- */
        PROCESS_MEMORY_COUNTERS info;
        GetProcessMemoryInfo( GetCurrentProcess( ), &info, sizeof(info) );
        if(screen == true)
                std::cout<<" MemoryPeak RSS:\t"<< (size_t)info.PeakWorkingSetSize <<" [B]\n \n";
        if(print == true)
                myfile << (size_t)info.PeakWorkingSetSize << "\n" ;
        return (size_t)info.PeakWorkingSetSize;

    #elif defined(__unix__) || defined(__unix) || defined(unix) || (defined(__APPLE__) && defined(__MACH__))
        /* BSD, Linux, and OSX -------------------------------------- */
        struct rusage rusage;
        getrusage( RUSAGE_SELF, &rusage );
        #if defined(__APPLE__) && defined(__MACH__)
            if(screen == true)
                std::cout<<" MemoryPeak RSS:\t"<< (size_t)rusage.ru_maxrss <<" [B]\n \n";
            if(print == true)
                myfile << (size_t)rusage.ru_maxrss << "\n" ;
            return (size_t)rusage.ru_maxrss;
        #else
            if(screen == true)
                std::cout<<" MemoryPeak RSS:\t"<< (size_t)(rusage.ru_maxrss * 1024L) <<" [B]\n \n";
            if(print == true)
                myfile << (size_t)(rusage.ru_maxrss * 1024L) << "\n" ;
            return (size_t)(rusage.ru_maxrss * 1024L);
        #endif

    #else
        /* Unknown OS ----------------------------------------------- */
        return (size_t)0L;			/* Unsupported. */

    #endif

    myfile.close();
}




/*
 * Returns the current resident set size (physical memory use) measured
 * in bytes, or zero if the value cannot be determined on this OS.
 */
size_t GetMemoryProcess(bool screen, bool print) 
{
    // open a file to write the memory consumption
    std::ofstream myfile;
    myfile.open ("MemoryProcess.txt", std::ofstream::out | std::ofstream::app);

    #if defined(_WIN32) || defined(WIN32) || defined(__WIN32) && !defined(__CYGWIN__)
        /* Windows -------------------------------------------------- */
        PROCESS_MEMORY_COUNTERS info;
        GetProcessMemoryInfo( GetCurrentProcess( ), &info, sizeof(info) );
        if(screen == true)
                std::cout<<" Memory RSS:\t"<< (size_t)info.WorkingSetSize <<" [B]\n \n";
        if(print == true)
                myfile << (size_t)info.WorkingSetSize << "\n" ;
        return (size_t)info.WorkingSetSize;

    #elif defined(__APPLE__) && defined(__MACH__)
        /* OSX ------------------------------------------------------ */
        struct mach_task_basic_info info;
        mach_msg_type_number_t infoCount = MACH_TASK_BASIC_INFO_COUNT;
        if ( task_info( mach_task_self( ), MACH_TASK_BASIC_INFO,
            (task_info_t)&info, &infoCount ) != KERN_SUCCESS )
            return (size_t)0L;		/* Can't access? */
        if(screen == true)
                std::cout<<" Memory RSS:\t"<< (size_t)info.resident_size <<" [B]\n \n";
        if(print == true)
                myfile << (size_t)info.resident_size << "\n" ;
        return (size_t)info.resident_size;

    #elif defined(__linux__) || defined(__linux) || defined(linux) || defined(__gnu_linux__)
        /* Linux ---------------------------------------------------- */
        long rss = 0L;
        FILE* fp = NULL;
        if ( (fp = fopen( "/proc/self/statm", "r" )) == NULL )
            return (size_t)0L;		/* Can't open? */
        if ( fscanf( fp, "%*s%ld", &rss ) != 1 )
        {
            fclose( fp );
            return (size_t)0L;		/* Can't read? */
        }
        fclose( fp );
        if(screen == true)
                std::cout<<" Memory RSS:\t"<< (size_t)rss * (size_t)sysconf( _SC_PAGESIZE) <<" [B]\n \n";
        if(print == true)
                myfile << (size_t)rss * (size_t)sysconf( _SC_PAGESIZE) << "\n" ;
        return (size_t)rss * (size_t)sysconf( _SC_PAGESIZE);

    #else
        /* AIX, BSD, Solaris, and Unknown OS ------------------------ */
        return (size_t)0L;			/* Unsupported. */
    #endif

    myfile.close();
}
