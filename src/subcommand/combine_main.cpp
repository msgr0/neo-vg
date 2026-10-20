/** \file combine_main.cpp
 *
 * Defines the "nvg combine" subcommand
 */


#include <omp.h>
#include <unistd.h>
#include <getopt.h>

#include <iostream>

#include <vg/io/vpkg.hpp>

#include "subcommand.hpp"

#include "../handle.hpp"
#include "../vg.hpp"
#include "../io/save_handle_graph.hpp"

using namespace std;
using namespace vg;
using namespace vg::subcommand;

void help_combine(char** argv) {
    cerr << "usage: " << argv[0] << " combine [options] <graph1.vg> [graph2.vg ...] >merged.vg" << endl
         << "Combines one or more graphs into a single file, regardless of input format." << endl
         << "Node IDs will be modified as needed to resolve conflicts (same as nvg ids -j)." << endl
         << endl
         << "Options:" << endl
         << "  -p, --connect-paths   add edges necessary to connect paths with the same name" << endl
         << "                        which are present in different graphs." << endl
         << "                        ex: If path x is present in graphs N-1 and N, then" << endl
         << "                        an edge connecting the last node of x in N-1 " << endl
         << "                        and the first node of x in N will be added." << endl
         << "  -h, --help            print this help message to stderr and exit" << endl;
}

int main_combine(int argc, char** argv) {
    Logger logger("nvg combine");

    if (argc == 2) {
        help_combine(argv);
        return 1;
    }

    bool connect_paths = false;
    
    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
        {
            {"help", no_argument, 0, 'h'},
            {"connect-paths", no_argument, 0, 'p'},
            {0, 0, 0, 0}
        };

        int option_index = 0;
        c = getopt_long (argc, argv, "h?p",
                         long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {
        case 'p':
            connect_paths = true;
            break;
        case 'h':
        case '?':
            help_combine(argv);
            exit(1);
            break;

        default:
            abort ();
        }
    }

    unique_ptr<MutablePathMutableHandleGraph> first_graph;
    string first_graph_filename = get_input_file_name(optind, argc, argv);
    first_graph = vg::io::VPKG::load_one<MutablePathMutableHandleGraph>(first_graph_filename);
    int64_t max_node_id = first_graph->max_node_id();

    while (optind < argc) {

        unique_ptr<MutablePathMutableHandleGraph> graph;
        string graph_filename = get_input_file_name(optind, argc, argv);
        graph = vg::io::VPKG::load_one<MutablePathMutableHandleGraph>(graph_filename);

        // join the id spaces if necessary
        int64_t delta = max_node_id - graph->min_node_id();
        if (delta >= 0) {
            graph->increment_node_ids(delta + 1);
        }
        max_node_id = graph->max_node_id();

        if (connect_paths) {
            handlealgs::append_path_handle_graph(graph.get(), first_graph.get(), true);
        } else {
            graph->for_each_path_handle([&](path_handle_t path_handle) {
                    string path_name = graph->get_path_name(path_handle);
                    if (first_graph->has_path(path_name)) {
                        logger.error() << "Paths with name \"" << path_name 
                                       << "\" found in multiple input graphs. If they are consecutive"
                                       << "subpath ranges, they can be connected by using the -p option." << endl;
                    }
                });
            handlealgs::copy_path_handle_graph(graph.get(), first_graph.get());
        }        
    }

    // Serialize the graph using VPKG.
    vg::io::save_handle_graph(first_graph.get(), cout);

    return 0;
}

// Register subcommand
static Subcommand vg_combine("combine", "merge multiple graph files together", main_combine);
