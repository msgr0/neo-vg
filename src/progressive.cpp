#include "progressive.hpp"

#include <vg/io/stream.hpp>
#include <iostream>
#include <iomanip>
#include <chrono>
#include <sys/ioctl.h>

namespace vg {

using namespace std;

void Progressive::with_progress(bool show_progress, const std::string& task, const std::function<void(const std::function<void(size_t, size_t)>& progress)>& callback) {
    if (!show_progress) {
        // Use the handy no-op function from libvgio.
        callback(vg::io::NO_PROGRESS);
    } else {
        // We really do need to show progress.
        Progressive progressive;
        progressive.show_progress = show_progress;
        progressive.preload_progress(task);
        bool first_progress_update = true;

        callback([&](size_t completed, size_t total) {
            if (completed != std::numeric_limits<size_t>::max() && total != std::numeric_limits<size_t>::max()) {
                // This is a real update;
                if (first_progress_update) {
                    // This is the first update. Make the bar.
                    progressive.create_progress(total);
                    first_progress_update = false;
                }
                // Tell the bar how big to be.
                progressive.update_progress(completed);
            }
        });

        progressive.destroy_progress();
    }
}

void Progressive::draw_progress_bar(const string& description, size_t done, size_t total, ostream& out) {
    // Work out how wide the terminal is.
    int console_width = 80;
    struct winsize win;
    if (ioctl(0, TIOCGWINSZ, &win) != -1) {
        console_width = win.ws_col;
    }

    if (done > total) {
        // Clear the line and complain.
        out << string(console_width, ' ') << "\r" << flush;
        cerr << "PROGRESS_BAR_EXCEPTION: _idx (" << done << ") went out of bounds, greater than n (" << total << ")." << endl << flush;
        return;
    }

    // Leave room for the description and the percentage.
    int bar_size = (int) ((console_width - (int) description.size() - 4) / 2.);
    double percent = total == 0 ? 100.0 : done * 100.0 / total;
    double percent_per_unit = 100.0 / bar_size;

    out << " " << description << " [";
    for (int i = 0; i < bar_size; i++) {
        out << (i * percent_per_unit < percent ? '=' : ' ');
    }
    out << "]" << setw(5) << setprecision(1) << fixed << percent << "%\r" << flush;
}

void Progressive::create_progress(const string& message, long count) {
    if (show_progress) {
        progress_message = message;
        create_progress(count);
    }
}

void Progressive::create_progress(long count) {
    if (show_progress) {
        progress_count = count;
        last_progress = 0;
        progress_seen = 0;
        progress_message.resize(30, ' ');
        progress = true;
        draw_progress_bar(progress_message, 0, progress_count);
    }
}

void Progressive::ensure_progress(long count) {
    if (show_progress) {
        if (!progress || progress_count != count) {
            create_progress(count);
        }
    }
}

void Progressive::preload_progress(const string& message) {
    if (show_progress && !progress) {
        // Set the message. We don't change it while a progress bar is
        // displayed.
        progress_message = message;
    }
}

void Progressive::update_progress(long i) {
    if (show_progress && progress) {
        if ((i <= progress_count
             && (long double) (i - last_progress) / (long double) progress_count >= 0.001)
            || i == progress_count) {
#pragma omp critical (progress)
            {
                draw_progress_bar(progress_message, i, progress_count);
                last_progress = i;
            }
        }
#pragma omp critical (progress_seen)
        {
            // Remember the amout of progress we have seen.
            progress_seen = i;
        }
    }
}

void Progressive::increment_progress() {
#pragma omp critical (increment_progress)
    {
        // Only one increment can happen at a time, so we don't get lower values
        // coming after higher values because the increments and the displays
        // interleaved.

        // What should we set the progress bar to?
        size_t new_progress;

#pragma omp critical (progress_seen)
        {   
            // We need to synchronize with update calls.
            
            // Bump up whatever we have to 1 more than that.
            new_progress = progress_seen + 1;
        }
        
        // TODO: interleaving updates and increments from different threads will
        // be weird, and progress may go backward. But that makes no real sense
        // anyway, so we just have to not explode/deadlock when it happens.
    
        update_progress(new_progress);
    }
}

void Progressive::destroy_progress(void) {
    if (show_progress && progress) {
        update_progress(progress_count);
        cerr << endl;
        progress_message = "progress";
        progress_count = 0;
    }
    progress = false;
}

}
