{
  description = "Build and test environment for neo-vg (nvg)";

  inputs.nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";

  outputs = { nixpkgs, ... }:
    let
      systems = [ "x86_64-linux" "aarch64-linux" "aarch64-darwin" "x86_64-darwin" ];
      forAllSystems = f: nixpkgs.lib.genAttrs systems (system: f nixpkgs.legacyPackages.${system});
    in {
      devShells = forAllSystems (pkgs:
        let
          # The Makefile looks for OpenMP (include/omp.h, lib/libomp) under one prefix
          openmp = pkgs.symlinkJoin {
            name = "openmp-prefix";
            paths = [ pkgs.llvmPackages.openmp pkgs.llvmPackages.openmp.dev ];
          };
        in {
          default = pkgs.mkShell {
            packages = with pkgs; [
              # Toolchain and build tools (also needed by the bundled deps' own builds)
              gnumake cmake pkg-config autoconf automake libtool gnum4 gettext
              bison flex gawk coreutils perl python3 git

              # Libraries nvg links against
              protobuf jansson boost openssl zlib bzip2 xz zstd lz4 ncurses
              gperftools python3Packages.pybind11

              # Tools used by the tests and the docs
              jq parallel bc curl wget unzip nodejs pandoc samtools bcftools htslib kmc
            ] ++ lib.optionals stdenv.isDarwin [ openmp ];

            # Absolute paths: CMake searches the Makefile's CMAKE_PREFIX_PATH (/usr) before PATH
            CC = "${pkgs.stdenv.cc}/bin/cc";
            CXX = "${pkgs.stdenv.cc}/bin/c++";

            # A system make found by CMake may be older than gnumake 4.4 and only
            # understand the pipe-style jobserver. nixpkgs has no static jansson
            # or abseil, so on Linux link system libraries dynamically (the
            # Makefile defaults to static there).
            MAKEFLAGS = "--jobserver-style=pipe"
              + pkgs.lib.optionalString pkgs.stdenv.isLinux " START_STATIC= END_STATIC=";
            HOMEBREW_PREFIX = pkgs.lib.optionalString pkgs.stdenv.isDarwin "${openmp}";
          };
        });
    };
}
