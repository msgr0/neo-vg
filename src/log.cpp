#include "log.hpp"
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <thread>
#include <vector>
#include <omp.h>

namespace vg {

/// Static initialization runs on the main thread.
static const std::thread::id main_thread_id = std::this_thread::get_id();

void fatal_exit() {
    if (omp_in_parallel() || std::this_thread::get_id() != main_thread_id) {
        // ponytail: skips atexit cleanup (temp files) too; a cooperative
        // cancel flag checked by the mapping loops would keep it.
        std::cerr.flush();
        fflush(stderr);
        std::_Exit(EXIT_FAILURE);
    }
    exit(EXIT_FAILURE);
}

cerrWrapper::~cerrWrapper() {
    std::string message = buffer.str();
    // One stdio call holds stderr's lock, so the message stays in one piece.
    std::cerr.flush();
    fwrite(message.data(), 1, message.size(), stderr);
    if (exit_on_destruct) {
        fatal_exit();
    }
}

namespace logging {

cerrWrapper info(const std::string& context) {
    return cerrWrapper("[" + context + "] ", false);
}

cerrWrapper warn(const std::string& context) {
    return cerrWrapper("warning[" + context + "] ", false);
}

cerrWrapper error(const std::string& context) {
    return cerrWrapper("error[" + context + "] ", true);
}

}

}
