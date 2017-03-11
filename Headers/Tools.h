///**************************************************************************
/// HEADER: Function Used As Tools (Memory, CPU,...)
///**************************************************************************

#ifndef TOOLS_H
#define TOOLS_H

#if defined(_WIN32) || defined(WIN32) || defined(__WIN32) && !defined(__CYGWIN__)
#include <windows.h>
#include <psapi.h>
    #pragma comment( lib, "psapi.lib" )

#elif defined(__unix__) || defined(__unix) || defined(unix) || (defined(__APPLE__) && defined(__MACH__))
#include <unistd.h>
#include <sys/resource.h>
    #include "sys/types.h"
    //#include "sys/sysinfo.h"

#if defined(__APPLE__) && defined(__MACH__)
#include <mach/mach.h>

#elif defined(__linux__) || defined(__linux) || defined(linux) || defined(__gnu_linux__)
#include <stdio.h>
    #include "sys/types.h"
    #include "sys/sysinfo.h"
    extern struct sysinfo memInfo;

#endif

#else
#error "Cannot define GetMemory( ) or GetMemoryProcessPeak( ) or GetMemoryProcess() for an unknown OS."
#endif

// Memory and CPU consumption
size_t GetMemory(bool screen, bool print);
size_t GetMemoryProcess(bool screen, bool print);
size_t GetMemoryProcessPeak(bool screen, bool print);

#endif
