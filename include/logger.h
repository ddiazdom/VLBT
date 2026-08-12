//
// Created by Diaz, Diego on 9.12.2024.
//
#ifndef VLBT_LOGGER_H
#define VLBT_LOGGER_H

#include <iostream>
#include <string>
#include <sstream>
#include <iomanip>
#include <chrono>
#ifdef USE_MALLOC_COUNT
#include <malloc_count.h>
#endif

//code generated with ChatGPT
#define CONCAT_IMPL(x,y) x##y
#define CONCAT(x,y) CONCAT_IMPL(x,y)

//from https://stackoverflow.com/questions/16605967/set-precision-of-stdto-string-when-converting-floating-point-values
template <typename T>
std::string to_string_with_precision(const T a_value, const int n = 6) {
    std::ostringstream out;
    out.precision(n);
    out << a_value;
    return std::move(out).str();
}

inline std::string format_time(double seconds) {
    std::ostringstream os;
    if (seconds < 60) {
        os << std::fixed << std::setprecision(3)<< seconds << "s";
    } else if (seconds < 3600) {
        int m = static_cast<int>(seconds) / 60;
        int s = static_cast<int>(seconds) % 60;
        os << m << "m"<< std::setw(2) << std::setfill('0') << s << "s";
    }
    else {
        int h = static_cast<int>(seconds) / 3600;
        int m = (static_cast<int>(seconds) % 3600) / 60;
        os << h << "h"<< std::setw(2) << std::setfill('0') << m << "m";
    }
    return os.str();
}

inline std::string format_space(off_t bytes){
    if(bytes<1000){
        return std::to_string(bytes)+" bytes";
    }

    if(bytes<1000000){
        float b = static_cast<float>(bytes)/1000;
        return to_string_with_precision(b, 3)+" KBs";
    }

    if(bytes < 1000000000L){
        float b = static_cast<float>(bytes)/1000000;
        return to_string_with_precision(b, 3)+" MBs";
    }

    if(bytes < 1000000000000L){
        double b = static_cast<double>(bytes)/1000000000L;
        return to_string_with_precision(b, 3)+" GBs";
    }

    double b = static_cast<double>(bytes)/1000000000000L;
    return to_string_with_precision(b, 3)+" TBs";
}

// -----------------------------
// Log levels
// -----------------------------
enum class log_level {
    ERROR = 0,
    WARN  = 1,
    INFO  = 2,
    DEBUG = 3,
    TRACE = 4
};

// -----------------------------
// Logger
// -----------------------------
struct Logger {

    static log_level level;
    static std::chrono::steady_clock::time_point start;
    static int indent;
    //static std::chrono::steady_clock::time_point last;

    static bool enabled(log_level msgLevel) {
        return static_cast<int>(msgLevel) <= static_cast<int>(level);
    }

    static void log(log_level msgLevel, const std::string& msg) {
        if (!enabled(msgLevel))
            return;

        std::ostream& out = (msgLevel<=log_level::WARN) ? std::cerr : std::cout;

        auto now = std::chrono::steady_clock::now();

        double total = std::chrono::duration<double>(now - start).count();
        //double delta = std::chrono::duration<double>(now - last).count();
        //last = now;

        int hours   = static_cast<int>(total) / 3600;
        int minutes = (static_cast<int>(total) % 3600) / 60;
        int seconds = static_cast<int>(total) % 60;

        out << "["
            << std::setw(2) << std::setfill('0') << hours << ":"
            << std::setw(2) << minutes << ":"
            << std::setw(2) << seconds << "] ";

        for (int i = 0; i < indent; ++i)
            out << "  ";
        out <<msg<<std::endl;
    }
};

int Logger::indent = 0;
log_level Logger::level = log_level::INFO;
std::chrono::steady_clock::time_point Logger::start = std::chrono::steady_clock::now();
//std::chrono::steady_clock::time_point Logger::last = start;

class ScopeIndent {
public:
    ScopeIndent() {
        Logger::indent++;
    }
    ~ScopeIndent() {
        Logger::indent--;
    }
};

// -----------------------------
// RAII Scope Logger
// -----------------------------
class ScopeLog {
public:
    ScopeLog(const std::string& name, log_level lvl)
        : name(name),
          level(lvl),
          active(Logger::enabled(lvl)),
          start(std::chrono::steady_clock::now()) {
        if (active) {
            Logger::log(level, "→ " + name);
#ifdef USE_MALLOC_COUNT
            malloc_count_reset_peak();
#endif
        }
    }

    ~ScopeLog() {
        if (!active) return;
        auto end = std::chrono::steady_clock::now();

        double elapsed = std::chrono::duration<double>(end - start).count();
        std::string msg = "← " + name + " (" + format_time(elapsed);
#ifdef USE_MALLOC_COUNT
        msg += ", " + format_space(static_cast<off_t>(malloc_count_peak()));
#endif
        msg += ")";
        Logger::log(level, msg);
    }

private:
    std::string name;
    log_level level;
    bool active;
    std::chrono::steady_clock::time_point start;
};

// -----------------------------
// Logging macros
// -----------------------------
#define TRACE_SCOPE() ScopeLog scope(__FUNCTION__, log_level::TRACE)
#define DEBUG_SCOPE() ScopeLog scope(__FUNCTION__, log_level::DEBUG)
#define SCOPE_INFO() ScopeIndent CONCAT(_scopeIndent_, __LINE__)

#define LOG_ERROR(msg) Logger::log(log_level::ERROR, msg)
#define LOG_WARN(msg)  Logger::log(log_level::WARN,  msg)
#define LOG_INFO(msg)  Logger::log(log_level::INFO,  msg)
#define LOG_DEBUG(msg) Logger::log(log_level::DEBUG, msg)
#define LOG_TRACE(msg) Logger::log(log_level::TRACE, msg)
#endif //VLBT_LOGGER_H
