#ifndef VG_LOG_HPP_INCLUDED
#define VG_LOG_HPP_INCLUDED

#include <string>
#include <iostream>
#include <sstream>

/** \file
 * log.hpp: defines a basic logging system for vg.
 * Shipped along with utility.hpp
 *
 * There are three options here:
 * - info(): do nothing but append a context string to the front of
 *   the log message. This is useful for logging progress messages.
 *   This is a generic logging utility.
 * - warn(): emit a warning message with a given context
 * - error(): emit an error message with a given context
 * 
 * Errors/warnings are output to std::cerr as:
 * error[<context>] message
 * warning[<context>] message
 * where "context" is a string like "nvg inject"
 * 
 * The error/warning functions exist to have a standardized error/warning format
 * across the vg subcommands. You don't *have* to use them elsewhere, but
 * it's a good idea to have consistent formatting.
 * 
 * These functions return a cerrWrapper object that behaves like a stream.
 * You can stream things to it with operator<<, and when the object goes out
 * of scope, it will write the message to std::cerr in one go and exit the program with a failure code (in the case of
 * error) or do nothing (in the case of warn() and info()).
 * 
 * Logging functions may either be accessed via logging::info/warn/error, or
 * via a Logger object that is initialized with a context string.
 */

namespace vg {

/// Exit with a failure code. Safe to call from worker threads: there it skips
/// static destructors, which would otherwise free data the other threads are
/// still using.
[[noreturn]] void fatal_exit();

// https://stackoverflow.com/a/25615354/
class cerrWrapper {
private:
    bool exit_on_destruct;
    // Buffered so the whole message reaches stderr in one write and does not
    // interleave with other threads' output.
    std::ostringstream buffer;

public:
    cerrWrapper(std::string prefix, bool exit_on_destruct) :
        exit_on_destruct(exit_on_destruct) {
        buffer << prefix;
    }

    template <typename T>
    cerrWrapper& operator<<(const T& t) {
        buffer << t;
        return *this;
    }

    cerrWrapper& operator<<(std::ostream& (*manip)(std::ostream&)) {
        buffer << manip;
        return *this;
    }

    ~cerrWrapper();
};

namespace logging {

/// Log to cerr with a standard format
/// "context" is caller context, e.g. "nvg inject"
cerrWrapper info(const std::string& context);
/// Emit a warning with a standard format
/// "context" is caller context, e.g. "nvg inject"
cerrWrapper warn(const std::string& context);
/// Error AND EXIT with a standard format
/// Once the cerrWrapper goes out of scope,
/// the program will exit with an error code.
/// "context" is caller context, e.g. "nvg inject"
cerrWrapper error(const std::string& context);

}

/// Class to set up at the start of a file
/// when you want to generate a bunch of loggers
class Logger {
private:
    std::string context;

public:
    Logger(const std::string& context) : context(context) {}

    /// Log to cerr with a standard format
    inline cerrWrapper info() const {
        return logging::info(context);
    }
    /// Emit a warning with a standard format
    inline cerrWrapper warn() const {
        return logging::warn(context);
    }
    /// Error AND EXIT with a standard format
    /// Once the cerrWrapper goes out of scope,
    /// the program will exit with an error code.
    inline cerrWrapper error() const {
        return logging::error(context);
    }
};
}

#endif
