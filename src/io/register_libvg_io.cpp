/**
 * \file register_libvg_io.hpp
 * Includes calls to register all libvg types with libvgio.
 */

#include "register_libvg_io.hpp"


namespace vg {

namespace io {

using namespace std;

// Each of these is defined in its own register_loader_*.cpp file.
void register_loader_saver_distance_index();
void register_loader_saver_gbwt();
void register_loader_saver_r_index();
void register_loader_saver_gbwtgraph();
void register_loader_saver_gbz();
void register_loader_saver_gbzgraph();
void register_loader_saver_gcsa();
void register_loader_saver_lcp();
void register_loader_saver_minimizer();
void register_loader_saver_snarl_manager();
void register_loader_saver_vg();
void register_loader_saver_xg();
void register_loader_saver_packed_graph();
void register_loader_saver_hash_graph();
void register_loader_saver_gfa();
void register_loader_saver_zip_codes();
void register_loader_gfaz();
void register_loader_params_json();

bool register_libvg_io() {
    register_loader_saver_distance_index();
    register_loader_saver_gbwt();
    register_loader_saver_r_index();
    register_loader_saver_gbwtgraph();
    register_loader_saver_gbz();
    register_loader_saver_gbzgraph();
    register_loader_saver_gcsa();
    register_loader_saver_lcp();
    register_loader_saver_minimizer();
    register_loader_saver_snarl_manager();
    register_loader_saver_vg();
    register_loader_saver_gfa();
    register_loader_gfaz();
    register_loader_saver_xg();
    register_loader_saver_packed_graph();
    register_loader_saver_hash_graph();
    register_loader_saver_zip_codes();
    register_loader_params_json();
    return true;
}
    
}

}
