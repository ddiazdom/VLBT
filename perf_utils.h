//
// Created by Diaz, Diego on 2.4.2025.
//

#ifndef BWT_DTS_BENCHMARKS_PERF_UTILS_H
#define BWT_DTS_BENCHMARKS_PERF_UTILS_H
#include <numeric>
#include <iostream>
#include <sstream>
#include <iomanip>
#include "performancecounters/benchmarker.h"

template<typename func>
int measure_avg_time(func function, std::vector<std::pair<uint8_t, uint64_t>>& tests, const uint16_t* vec, std::string name) {

    std::vector<long long> durations;

    for (auto & test : tests) {
        auto start = std::chrono::high_resolution_clock::now();
        function(vec, test.first, test.second);
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
        durations.push_back(duration);
    }

    // Compute average duration
    long long totalDuration = std::accumulate(durations.begin(), durations.end(), 0LL);
    double averageDuration = static_cast<double>(totalDuration) / static_cast<double>(tests.size());
    std::cout << "Average running time "<<name<<": "<< averageDuration << " nanoseconds" << std::endl;
    return 0;
}

void pretty_print(size_t volume, size_t bytes, std::string name,
                  event_aggregate agg) {
    printf("%-40s : ", name.c_str());
    printf(" %5.2f GB/s ", bytes / agg.fastest_elapsed_ns());
    printf(" %5.1f Ma/s ", volume * 1000.0 / agg.fastest_elapsed_ns());
    printf(" %5.2f ns/d ", agg.fastest_elapsed_ns() / volume);
    if (collector.has_events()) {
        printf(" %5.2f GHz ", agg.fastest_cycles() / agg.fastest_elapsed_ns());
        printf(" %5.2f c/d ", agg.fastest_cycles() / volume);
        printf(" %5.2f i/d ", agg.fastest_instructions() / volume);
        printf(" %5.2f c/b ", agg.fastest_cycles() / bytes);
        printf(" %5.2f i/b ", agg.fastest_instructions() / bytes);
        printf(" %5.2f i/c ", agg.fastest_instructions() / agg.fastest_cycles());
        printf(" %5.2f br_miss ", agg.branch_misses());
    }
    printf("\n");
}

template<class time_t>
std::string report_time(time_t start, time_t end, size_t padding){
    auto dur = end - start;
    auto h = std::chrono::duration_cast<std::chrono::hours>(dur);
    auto m = std::chrono::duration_cast<std::chrono::minutes>(dur -= h);
    auto s = std::chrono::duration_cast<std::chrono::seconds>(dur -= m);
    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(dur -= s);
    std::stringstream time;

    for(size_t i=0;i<padding;i++) std::cout<<" ";
    if(h.count()>0){
        time <<"build_time (hh:mm:ss.ms): "<<std::setfill('0')<<std::setw(2)<<h.count()<<":"<<std::setfill('0')<<std::setw(2)<<m.count()<<":"<<std::setfill('0')<<std::setw(2)<<s.count()<<"."<<ms.count();
    }else if(m.count()>0){
        time <<"build_time (mm:ss.ms): "<<std::setfill('0')<<std::setw(2)<<size_t(m.count())<<":"<<std::setfill('0')<<std::setw(2)<<s.count()<<"."<<ms.count();
    }else if(s.count()>0){
        time <<"build_time (ss.ms): "<<std::setfill('0')<<std::setw(2)<<size_t(s.count())<<"."<<ms.count();
    }else{
        time <<"build_time (ms): "<<ms.count();
    }

    return time.str();
}
#endif //BWT_DTS_BENCHMARKS_PERF_UTILS_H
