/** \file help_main.cpp
 *
 * Defines the "nvg help" subcommand, which describes subcommands.
 */


#include <omp.h>
#include <unistd.h>
#include <getopt.h>
#include <stdlib.h>

#include <iostream>
#include <string>
#include <vector>

#include "subcommand.hpp"
#include "../version.hpp"

using namespace std;
using namespace vg;
using namespace vg::subcommand;

namespace vg {
namespace subcommand {

/// Print the version banner, then every subcommand in each of the given
/// categories (names padded to name_width), then the footer text and the
/// support line. Also used by main.cpp when nvg is run with no command.
void print_subcommand_help(char** argv, const vector<CommandCategory>& categories,
                           size_t name_width, const string& footer) {
    cerr << "nvg: variation graph tool, version " << Version::get_short() << endl
         << endl
         << "usage: " << argv[0] << " <command> [options]" << endl
         << endl;

    for (auto category : categories) {
        cerr << category << ":" << endl;
        Subcommand::for_each(category, [&](const Subcommand& command) {
            // Pad all the names so the descriptions line up
            string name = command.get_name();
            name.resize(name_width, ' ');
            cerr << "  -- " << name << command.get_description() << endl;
        });
        cerr << endl;
    }

    cerr << footer
         << "For technical support, please visit: https://www.biostars.org/tag/vg/" << endl << endl;
}

}
}

int main_help(int argc, char** argv) {
    print_subcommand_help(argv, {PIPELINE, TOOLKIT, WIDGET, DEVELOPMENT}, 18, "");
    return 0;
}

// Register subcommand
static Subcommand vg_help("help", "show all subcommands", PIPELINE, main_help);

