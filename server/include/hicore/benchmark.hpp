// time recording from https://stackoverflow.com/questions/17432502/how-can-i-measure-cpu-time-and-wall-clock-time-on-both-linux-windows
#ifndef BENCHMARK_HPP_
#define BENCHMARK_HPP_
#include <iostream>
//  Windows
#ifdef _WIN32
#include <Windows.h>
inline double get_wall_time(){
    LARGE_INTEGER time,freq;
    if (!QueryPerformanceFrequency(&freq)){
        //  Handle error
        return 0;
    }
    if (!QueryPerformanceCounter(&time)){
        //  Handle error
        return 0;
    }
    return (double)time.QuadPart / freq.QuadPart;
}
inline double get_cpu_time(){
    FILETIME a,b,c,d;
    if (GetProcessTimes(GetCurrentProcess(),&a,&b,&c,&d) != 0){
        //  Returns total user time.
        //  Can be tweaked to include kernel times as well.
        return
            (double)(d.dwLowDateTime |
            ((unsigned long long)d.dwHighDateTime << 32)) * 0.0000001;
    }else{
        //  Handle error
        return 0;
    }
}

//  Posix/Linux
#else
#include <time.h>
#include <sys/time.h>
inline double get_wall_time(){
    struct timeval time;
    if (gettimeofday(&time,NULL)){
        //  Handle error
        return 0;
    }
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
}
inline double get_cpu_time(){
    return (double)clock() / CLOCKS_PER_SEC;
}
#endif

#define BENCH_BEGIN(argc, argv)                                                                    \
    for (int i = 0; i < argc; i++) {                                                               \
        if (i > 0)                                                                                 \
            std::cerr << " ";                                                                      \
        std::cerr << argv[i];                                                                      \
    }                                                                                              \
    std::cerr << ",";                                                                              \
    auto __bench_t__ = get_wall_time();                                                            \
    std::vector<double> __bench_ts__;
#define BENCH_SECTION                                                                              \
    if (!__bench_ts__.empty()) {                                                                   \
        __bench_ts__[__bench_ts__.size() - 1] = get_wall_time() - __bench_ts__.back();             \
    }                                                                                              \
    __bench_ts__.push_back(get_wall_time());
#define BENCH_END                                                                                  \
    if (!__bench_ts__.empty()) {                                                                   \
        __bench_ts__[__bench_ts__.size() - 1] = get_wall_time() - __bench_ts__.back();             \
    }                                                                                              \
    for (auto &t : __bench_ts__) {                                                                 \
        std::cerr << std::fixed << t << ",";                                                       \
    }                                                                                              \
    std::cerr << std::fixed << (get_wall_time() - __bench_t__) << std::endl;

#endif // BENCHMARK_HPP_