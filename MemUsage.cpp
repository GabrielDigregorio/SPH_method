#include "SPH.hpp"

#include <windows.h>
#include <stdio.h>
#include <psapi.h>

// For Windows
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__)
    
    #include "windows.h"
    #include <stdio.h>
    #include "psapi.h"

    // Get the memory usage of the current process
    void GetMemoryProcess(bool print){

    }

    // Get the memory usage of the computer
    void GetMemory(bool print){

        // open a file to write the memory consumption
        std::ofstream myfile;
        myfile.open ("MemoryConsumption.txt", std::ofstream::out | std::ofstream::app);

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

        if(print == true){
            std::cout<<"\n TotVirtMem:    \t"<< totalVirtualMem <<" [KB]\n";
            std::cout<<" VirtMemCurrUsed:\t"<< virtualMemUsed <<" [KB]\n";
            std::cout<<" TotPhysMem(RAM):\t"<< totalPhysMem <<" [KB]\n";
            std::cout<<" PhysMemCurrUsed:\t"<< physMemUsed <<" [KB]\n \n";
        }
            
        myfile << totalVirtualMem << " " << virtualMemUsed << " " << totalPhysMem << " " << physMemUsed << "\n" ;
        myfile.close();

    }

// For Linux & UNIX
#elif __unix || __unix__ || __linux__
    #include "stdlib.h"
    #include "stdio.h"
    #include "string.h"
    #include "sys/types.h"
    #include "sys/sysinfo.h"

    struct sysinfo memInfo;

    // Get the memory usage of the current process
    void GetMemoryProcess(bool print){

    }

    // Get the memory usage of the computer
    void GetMemory(bool print){

        // open a file to write the memory consumption
        std::ofstream myfile;
        myfile.open ("MemoryConsumption.txt", std::ofstream::out | std::ofstream::app);

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

        if(print == true){
            std::cout<<"\n TotVirtMem:    \t"<< totalVirtualMem <<" [KB]\n";
            std::cout<<" VirtMemCurrUsed:\t"<< virtualMemUsed <<" [KB]\n";
            std::cout<<" TotPhysMem(RAM):\t"<< totalPhysMem <<" [KB]\n";
            std::cout<<" PhysMemCurrUsed:\t"<< physMemUsed <<" [KB]\n \n";
        }
            
        myfile << totalVirtualMem << " " << virtualMemUsed << " " << totalPhysMem << " " << physMemUsed << "\n" ;
        myfile.close();
    }

    int parseLine(char* line){
        // This assumes that a digit will be found and the line ends in " Kb".
        int i = strlen(line);
        const char* p = line;
        while (*p <'0' || *p > '9') p++;
        line[i-3] = '\0';
        i = atoi(p);
        return i;
    }

    //Note: this value is in KB!
    int getValue(){
        FILE* file = fopen("/proc/self/status", "r");
        int result = -1;
        char line[128];

        while (fgets(line, 128, file) != NULL){
            if (strncmp(line, "VmSize:", 7) == 0){
                result = parseLine(line);
                break;
            }
        }
        fclose(file);
        return result;
    }

#endif
