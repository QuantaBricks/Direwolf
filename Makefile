# Copyright (c) 2026 QuantaBricks
# SPDX-License-Identifier: Apache-2.0

# Checked at Makefile-parse time (not as a recipe/prerequisite) so it's
# immune to `make -j` scheduling and always fails fast, before wasting
# time partway into a third-party cmake build with a cryptic "cmake:
# not found" / Error 127.
ifeq ($(shell command -v cmake 2>/dev/null),)
$(error cmake not found - install it first (e.g. 'sudo apt install cmake' on Debian/Ubuntu), then re-run make)
endif
ifeq ($(shell command -v gfortran 2>/dev/null),)
$(error gfortran not found - install it first (e.g. 'sudo apt install gfortran' on Debian/Ubuntu), then re-run make)
endif

# Override for a portable release binary, e.g. `make MARCH=x86-64-v2`:
# a binary built with the default -march=native only runs on CPUs with
# the exact same feature set as the build machine (illegal instruction
# otherwise) - fine for a recipient building their own copy, not fine
# for a prebuilt binary handed to someone else's machine.
MARCH      ?= native

FORT90     = gfortran
CC         = gcc
DATADIR    = data
LIBDIR     = lib
BUILDDIR   = build
EXDIR      = examples
NPROC      = $(shell nproc)
# third_party/* is all vendored source, built automatically by `make`
# (never a prebuilt binary) - see README.md's "Third-party dependencies"
# table for licenses and third_party/OpenBLAS's note on why OpenBLAS is
# built USE_THREAD=1 (thread-safety, not just convenience).
LIBCINTDIR    = third_party/libcint
LIBCINTECPDIR = third_party/libcint_ecp
LIBCINTINCDIR = $(LIBCINTDIR)/build/include
tagCOMPC   = -c -O3 -I$(LIBCINTINCDIR) -I$(LIBCINTECPDIR)
DFTD3DIR   = third_party/simple-dftd3
DFTD4DIR   = third_party/dftd4
MULTICHARGEDIR = third_party/multicharge
GCPDIR     = third_party/gcp
MCTCDIR    = third_party/mctc-lib
TOMLFDIR   = third_party/toml-f
THIRDPARTY_INSTALL = third_party/_install
tagCOMPD3  = -c -free -O3 -fopenmp -I$(DFTD3DIR)/build/include -I$(MCTCDIR)/build/include -J$(BUILDDIR)
tagCOMPD4  = -c -free -O3 -fopenmp -I$(DFTD4DIR)/build/include -I$(MCTCDIR)/build/include -J$(BUILDDIR)
tagCOMPGCP = -c -free -O3 -fopenmp -I$(GCPDIR)/build/include -J$(BUILDDIR)
OPENBLASDIR = third_party/OpenBLAS
LIBXCDIR   = third_party/libxc
tagCOMP    = -c -free -O3 -fopenmp -I$(DATADIR) -I$(LIBXCDIR)/build -J$(BUILDDIR)
tagCOMPF77 = -c -O3 -fopenmp -I$(DATADIR) -J$(BUILDDIR)

VPATH = src/main:src/starter:src/starter/io:src/starter/io/namelist:src/starter/io/xyz:src/starter/io/column:src/core:src/integrals:src/integrals/core:src/integrals/exact:src/integrals/df:src/integrals/cosx:src/integrals/force:src/integrals/hessian:src/localization:src/xc:src/dispersion:src/scf:src/force:src/guess:src/util:src/solvent:src/solvent/cavity:src/solvent/smd:src/solvent/cosmo_impl:examples

objects = $(BUILDDIR)/mod_data.o \
          $(BUILDDIR)/mod_profile.o \
          $(BUILDDIR)/mod_scf_history.o \
          $(BUILDDIR)/mod_meminfo.o \
          $(BUILDDIR)/mod_mem_predict.o \
          $(BUILDDIR)/exe_path_shim.o \
          $(BUILDDIR)/mod_basis_files.o \
          $(BUILDDIR)/mod_integrals.o \
          $(BUILDDIR)/mod_integrals_core.o \
          $(BUILDDIR)/mod_integrals_atomic_guess.o \
          $(BUILDDIR)/mod_integrals_puream.o \
          $(BUILDDIR)/properties.o \
          $(BUILDDIR)/mod_integrals_coulomb.o \
          $(BUILDDIR)/mod_integrals_exchange.o \
          $(BUILDDIR)/mod_integrals_exchange_cosx.o \
          $(BUILDDIR)/mod_integrals_exchange_cosx_k.o \
          $(BUILDDIR)/mod_integrals_exchange_cosx_lr.o \
          $(BUILDDIR)/mod_integrals_exchange_cosx_sr.o \
          $(BUILDDIR)/mod_integrals_exchange_cosx_force.o \
          $(BUILDDIR)/mod_integrals_force.o \
          $(BUILDDIR)/mod_integrals_hessian.o \
          $(BUILDDIR)/hessian_1e.o \
          $(BUILDDIR)/mod_integrals_df_setup.o \
          $(BUILDDIR)/mod_integrals_df_coulomb.o \
          $(BUILDDIR)/mod_integrals_df_exchange.o \
          $(BUILDDIR)/mod_integrals_df_force_store.o \
          $(BUILDDIR)/mod_integrals_df_force_direct_j.o \
          $(BUILDDIR)/mod_integrals_df_force_direct_jk.o \
          $(BUILDDIR)/mod_exchange.o \
          $(BUILDDIR)/mod_density_fitting.o \
          $(BUILDDIR)/mod_localization.o \
          $(BUILDDIR)/mod_xc.o \
          $(BUILDDIR)/mod_vv10.o \
          $(BUILDDIR)/mod_vv10_aca.o \
          $(BUILDDIR)/mod_dispersion.o \
          $(BUILDDIR)/mod_dispersion_d3.o \
          $(BUILDDIR)/mod_dispersion_d4.o \
          $(BUILDDIR)/mod_gcp.o \
          $(BUILDDIR)/mod_cosmo_constants.o \
          $(BUILDDIR)/mod_cosmo_dual.o \
          $(BUILDDIR)/mod_cosmo_radii.o \
          $(BUILDDIR)/mod_cosmo_geodesic.o \
          $(BUILDDIR)/mod_cosmo_clip.o \
          $(BUILDDIR)/mod_cosmo_gaubon.o \
          $(BUILDDIR)/mod_cosmo_cavity.o \
          $(BUILDDIR)/mod_smd_tables.o \
          $(BUILDDIR)/mod_smd_sts.o \
          $(BUILDDIR)/mod_smd_daareal.o \
          $(BUILDDIR)/mod_smd_cds.o \
          $(BUILDDIR)/mod_integrals_cosmo.o \
          $(BUILDDIR)/mod_cosmo_state.o \
          $(BUILDDIR)/mod_cosmo_solvents.o \
          $(BUILDDIR)/mod_cosmo_init.o \
          $(BUILDDIR)/mod_cosmo_scf.o \
          $(BUILDDIR)/mod_cosmo_force.o \
          $(BUILDDIR)/cosmo.o \
          $(BUILDDIR)/tool_io.o \
          $(BUILDDIR)/mod_checkpoint.o \
          $(BUILDDIR)/mod_molden.o \
          $(BUILDDIR)/guess.o \
          $(BUILDDIR)/SCF.o \
          $(BUILDDIR)/scf_adiis.o \
          $(BUILDDIR)/scf_diis.o \
          $(BUILDDIR)/scf_fock.o \
          $(BUILDDIR)/solveKS.o \
          $(BUILDDIR)/cphf.o \
          $(BUILDDIR)/tool_mat.o \
          $(BUILDDIR)/grid_gen.o \
          $(BUILDDIR)/grid_dynamic.o \
          $(BUILDDIR)/gto_eval.o \
          $(BUILDDIR)/Lebedev.o \
          $(BUILDDIR)/DFT.o \
          $(BUILDDIR)/DFT_vv10_pass.o \
          $(BUILDDIR)/DFT_shellpair.o \
          $(BUILDDIR)/force.o \
          $(BUILDDIR)/engine.o

# File-driven-frontend-only IO (src/starter/io): namelist/&atoms-file parsing
# for the example/regression driver(s) below (run_engine, and any
# future run_opt/run_bomd). Deliberately NOT part of $(objects)/
# libengine.a - that archive is the in-process integration API
# (EngineUp itself, called directly e.g. via iso_c_binding from
# Python), which never touches a file on disk, so it has no business
# depending on a file-format parser. See mod_engine_input.f90's header.
io_objects = $(BUILDDIR)/mod_engine_input_types.o \
             $(BUILDDIR)/mod_engine_input_elements.o \
             $(BUILDDIR)/mod_engine_input_block.o \
             $(BUILDDIR)/mod_engine_input_atoms_namelist.o \
             $(BUILDDIR)/mod_engine_input_atoms_xyz.o \
             $(BUILDDIR)/mod_engine_input_atoms_column.o \
             $(BUILDDIR)/mod_engine_input.o


all: Engine

$(BUILDDIR):
	mkdir -p $(BUILDDIR)

# Pre-created eagerly so every -I$(LIBXCDIR)/build compile flag resolves to a
# real directory even before libxc's own cmake build has run (order-only
# prereq below) - otherwise a parallel `make -j` compiles unrelated .f90
# files before that directory exists and gfortran spams "Nonexistent
# include directory" warnings. cmake -B build (libxc's own rule) is fine
# with the directory already existing.
$(LIBXCDIR)/build:
	mkdir -p $(LIBXCDIR)/build

# mod_version.f90 is regenerated on every build (.PHONY, never up-to-date)
# so the embedded git hash always reflects the tree actually being built,
# not whatever commit happened to be HEAD the last time this file was
# generated. Cheap (one git call), so unconditional regeneration is fine.
# The major.minor release number lives in VERSION (repo root) - bump it
# there once per release, not here.
$(BUILDDIR)/mod_version.f90: VERSION | $(BUILDDIR)
	@git_hash=$$(git rev-parse --short=12 HEAD 2>/dev/null || echo unknown); \
	engine_ver=$$(cat VERSION); \
	printf 'module mod_version\nimplicit none\ncharacter(len=*),parameter :: engine_version = "%s"\ncharacter(len=*),parameter :: engine_git_version = "%s"\nend module mod_version\n' "$$engine_ver" "$$git_hash" > $@
.PHONY: $(BUILDDIR)/mod_version.f90

$(BUILDDIR)/mod_version.o: $(BUILDDIR)/mod_version.f90
	$(FORT90) $(tagCOMP) $< -o $@

# VV10's O(npts^2) double sum is division-throughput-bound; -funsafe-math-
# optimizations lets the 3 reciprocal-per-pair become fast approximate rcp
# instructions (~2 ulp), a controlled numerical trade for ~2x speedup in the
# dominant NLC pass (verified against exact division: NLC energy agrees to
# <1e-11 Ha). mod_vv10_aca's phi4 kernel uses the same real(4)/single-
# division form, so it gets the same flags for the same reason.
$(BUILDDIR)/mod_vv10_aca.o: src/xc/mod_vv10_aca.f90 | $(BUILDDIR)
	$(FORT90) $(tagCOMP) -funsafe-math-optimizations -march=$(MARCH) $< -o $@

$(BUILDDIR)/mod_vv10.o: src/xc/mod_vv10.f90 $(BUILDDIR)/mod_data.o | $(BUILDDIR)
	$(FORT90) $(tagCOMP) -funsafe-math-optimizations -march=$(MARCH) $< -o $@

$(BUILDDIR)/mod_dispersion_d3.o: src/dispersion/mod_dispersion_d3.f90 $(DFTD3DIR)/build/libs-dftd3.so | $(BUILDDIR)
	$(FORT90) $(tagCOMPD3) $< -o $@

$(BUILDDIR)/mod_dispersion_d4.o: src/dispersion/mod_dispersion_d4.f90 $(DFTD4DIR)/build/libdftd4.so | $(BUILDDIR)
	$(FORT90) $(tagCOMPD4) $< -o $@

$(BUILDDIR)/mod_gcp.o: src/dispersion/mod_gcp.f90 $(GCPDIR)/build/libgcp.so | $(BUILDDIR)
	$(FORT90) $(tagCOMPGCP) $< -o $@

$(BUILDDIR)/%.o: %.f90 | $(BUILDDIR)
	$(FORT90) $(tagCOMP) $< -o $@

$(BUILDDIR)/%.o: %.F | $(BUILDDIR)
	$(FORT90) $(tagCOMPF77) $< -o $@

# --- third-party dependency chain: built automatically on first `make` ---
$(LIBCINTDIR)/build/libcint.a:
	cd $(LIBCINTDIR) && cmake -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=0 && cmake --build build -j$(NPROC)

$(OPENBLASDIR)/libopenblas.a:
	$(MAKE) -C $(OPENBLASDIR) USE_THREAD=1 USE_OPENMP=0 NO_SHARED=1 -j$(NPROC)

# -DBUILD_TESTING=OFF skips libxc's own testsuite/ subdirectory (the
# xc-regression/xc-consistency/xc-threshold/xc-info/xc-get_data
# executables + decompressing/diffing every .bz2 regression case under
# third_party/libxc/testsuite/) - none of that is needed to produce
# libxcf03.a/libxc.a, it just makes every from-scratch build slower.
$(LIBXCDIR)/build/libxcf03.a:
	cd $(LIBXCDIR) && cmake -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=0 -DENABLE_FORTRAN=ON -DENABLE_XHOST=OFF -DBUILD_TESTING=OFF && cmake --build build -j$(NPROC)

# -DBUILD_TESTING=OFF: toml-f's CMakeLists (via `include(CTest)`) tries
# to FetchContent the test-drive framework from GitHub when its own
# test suite is enabled (default ON) - unneeded to produce
# libtoml-f.a, and a real failure mode: on a machine where that
# FetchContent step doesn't cleanly complete, cmake configure fails
# outright ("add_subdirectory given source 'build' which is not an
# existing directory") and takes the whole build down with it.
$(TOMLFDIR)/build/libtoml-f.a:
	cd $(TOMLFDIR) && cmake -B build -DCMAKE_PREFIX_PATH=$(abspath $(THIRDPARTY_INSTALL)) -DCMAKE_INSTALL_PREFIX=$(abspath $(THIRDPARTY_INSTALL)) -DBUILD_TESTING=OFF && cmake --build build -j$(NPROC) && cmake --install build

$(MCTCDIR)/build/libmctc-lib.a: $(TOMLFDIR)/build/libtoml-f.a
	cd $(MCTCDIR) && cmake -B build -DCMAKE_PREFIX_PATH=$(abspath $(THIRDPARTY_INSTALL)) -DCMAKE_INSTALL_PREFIX=$(abspath $(THIRDPARTY_INSTALL)) && cmake --build build -j$(NPROC) && cmake --install build

$(DFTD3DIR)/build/libs-dftd3.so: $(MCTCDIR)/build/libmctc-lib.a
	cd $(DFTD3DIR) && cmake -B build -DCMAKE_PREFIX_PATH=$(abspath $(THIRDPARTY_INSTALL)) -DCMAKE_INSTALL_PREFIX=$(abspath $(THIRDPARTY_INSTALL)) -DBUILD_SHARED_LIBS=ON && cmake --build build -j$(NPROC)

$(GCPDIR)/build/libgcp.so:
	cd $(GCPDIR) && cmake -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=ON && cmake --build build -j$(NPROC)

# -DBLA_VENDOR=OpenBLAS + -DCMAKE_LIBRARY_PATH=<vendored OpenBLAS dir>:
# multicharge's CMakeLists independently requires a LAPACK/BLAS
# (find_package("custom-lapack" REQUIRED), unrelated to the
# already-vendored-and-built $(OPENBLASDIR) every other target here
# links against) - without this it demands a SYSTEM libblas-dev/
# liblapack-dev be installed, defeating the entire point of vendoring
# OpenBLAS from source (a from-scratch build on a machine with no
# system BLAS/LAPACK installed fails at this step otherwise - found on
# a real from-scratch test on a machine that had neither).
$(MULTICHARGEDIR)/build/libmulticharge.a: $(MCTCDIR)/build/libmctc-lib.a $(OPENBLASDIR)/libopenblas.a
	cd $(MULTICHARGEDIR) && cmake -B build -DCMAKE_PREFIX_PATH=$(abspath $(THIRDPARTY_INSTALL)) -DCMAKE_INSTALL_PREFIX=$(abspath $(THIRDPARTY_INSTALL)) -DBLA_VENDOR=OpenBLAS -DCMAKE_LIBRARY_PATH=$(abspath $(OPENBLASDIR)) && cmake --build build -j$(NPROC) && cmake --install build

# Same independent-BLAS issue as multicharge above (dftd4's own
# find_package("custom-blas") REQUIRED) - same fix.
$(DFTD4DIR)/build/libdftd4.so: $(MULTICHARGEDIR)/build/libmulticharge.a $(OPENBLASDIR)/libopenblas.a
	cd $(DFTD4DIR) && cmake -B build -DCMAKE_PREFIX_PATH=$(abspath $(THIRDPARTY_INSTALL)) -DCMAKE_INSTALL_PREFIX=$(abspath $(THIRDPARTY_INSTALL)) -DBUILD_SHARED_LIBS=ON -DBLA_VENDOR=OpenBLAS -DCMAKE_LIBRARY_PATH=$(abspath $(OPENBLASDIR)) && cmake --build build -j$(NPROC)
# --- end third-party dependency chain ---

$(BUILDDIR)/nr_ecp.o: $(LIBCINTECPDIR)/nr_ecp.c $(LIBCINTDIR)/build/libcint.a | $(BUILDDIR)
	$(CC) $(tagCOMPC) $< -o $@

$(BUILDDIR)/nr_ecp_deriv.o: $(LIBCINTECPDIR)/nr_ecp_deriv.c $(LIBCINTDIR)/build/libcint.a | $(BUILDDIR)
	$(CC) $(tagCOMPC) $< -o $@

$(BUILDDIR)/ecp2f.o: $(LIBCINTECPDIR)/ecp2f.c $(LIBCINTDIR)/build/libcint.a | $(BUILDDIR)
	$(CC) $(tagCOMPC) $< -o $@

# exe_path_shim.c: readlink("/proc/self/exe",...) wrapper, called from
# mod_basis_files.f90 so data/bases/ resolves relative to the running
# binary's own location, not the caller's cwd - see that file's header.
$(BUILDDIR)/exe_path_shim.o: src/core/exe_path_shim.c | $(BUILDDIR)
	$(CC) -c -O3 $< -o $@

$(LIBDIR):
	mkdir -p $(LIBDIR)

$(LIBDIR)/libcint_ecp.a: $(BUILDDIR)/nr_ecp.o $(BUILDDIR)/nr_ecp_deriv.o $(BUILDDIR)/ecp2f.o | $(LIBDIR)
	ar rcs $@ $^

# module-dependency ordering (explicit prereqs so .mod files exist before use)
$(BUILDDIR)/mod_profile.o   : $(BUILDDIR)/mod_meminfo.o
$(BUILDDIR)/mod_scf_history.o : $(BUILDDIR)/mod_meminfo.o
$(BUILDDIR)/mod_integrals.o : $(BUILDDIR)/mod_data.o $(BUILDDIR)/mod_meminfo.o $(BUILDDIR)/mod_profile.o
# mod_integrals.f90 declares the mod_integrals module's public interface
# only (a separate-module-procedure "interface ... end interface" block,
# no bodies); the files below are Fortran SUBMODULES of it (each
# "submodule (mod_integrals) <name>_impl") holding the actual
# implementations, split by feature purely for file-size/readability
# reasons - libcint calls still only ever appear inside this module
# family, matching mod_integrals.f90's own "sole boundary" header
# comment. Each needs mod_integrals.o built first (for its .mod file);
# no ordering is needed AMONG them, they don't depend on each other's
# .mod files (only on mod_integrals.mod, host-associated).
#
# src/integrals/df/ (2026-07-25 reorganization) is the RI-J/RI-K density-
# fitting family specifically - was one 3347-line mod_integrals_df.f90,
# split by feature the same way core/coulomb/exchange/force already were:
# setup (aux basis + compact-pair bookkeeping), coulomb (RI-J STORE+
# DIRECT), exchange (RI-K STORE+DIRECT), force_store (STORE-mode blocked
# JK gradient path), force_direct_j/force_direct_jk (DIRECT-mode gradient
# paths). df_force_2c2e_contract is a module subroutine (declared in
# mod_integrals.f90) rather than living in just one of these files since
# all three force files call it.
#
# 2026-08-13 follow-up reorganization: the same mod_integrals.f90 submodule
# family had grown a second axis (COSX semi-numerical K alongside exact
# STORE/DIRECT and DF/RI-J/RI-K, plus force/gradient and Hessian code for
# every one of those build modes) all sitting flat in src/integrals/ and
# src/force/ - pure file moves (git mv, no code changed) into subfolders
# by that axis: src/integrals/core/ (init/Schwarz/shell-pairs + the
# libcint kernel wrapper layer), src/integrals/exact/ (STORE+DIRECT J/K -
# each file is still its own STORE+DIRECT dispatcher, see its own header;
# NOT split further per-mode since STORE/DIRECT share pk_list/screening
# helpers within a file), src/integrals/cosx/, src/integrals/force/
# (gradients for exact + all of DF, moved out of src/integrals/ root and
# src/integrals/df/), src/integrals/hessian/ (moved out of src/force/,
# which was always just a stray placement - these are submodule
# (mod_integrals) files, not part of force.f90's own driver code).
# src/force/force.f90 itself stays put: it's the SCF-level force
# orchestration/assembly logic, not a mod_integrals implementation file.
$(BUILDDIR)/mod_integrals_core.o     : $(BUILDDIR)/mod_integrals.o $(BUILDDIR)/mod_data.o $(BUILDDIR)/mod_meminfo.o $(BUILDDIR)/mod_mem_predict.o $(BUILDDIR)/mod_exchange.o $(BUILDDIR)/mod_basis_files.o
$(BUILDDIR)/mod_integrals_atomic_guess.o : $(BUILDDIR)/mod_integrals.o $(BUILDDIR)/mod_data.o
$(BUILDDIR)/mod_integrals_puream.o : $(BUILDDIR)/mod_integrals.o
$(BUILDDIR)/properties.o : $(BUILDDIR)/mod_integrals.o $(BUILDDIR)/mod_data.o
$(BUILDDIR)/mod_integrals_coulomb.o  : $(BUILDDIR)/mod_integrals.o
$(BUILDDIR)/mod_integrals_exchange.o : $(BUILDDIR)/mod_integrals.o $(BUILDDIR)/mod_exchange.o
$(BUILDDIR)/mod_integrals_exchange_cosx.o : $(BUILDDIR)/mod_integrals.o $(BUILDDIR)/mod_exchange.o $(BUILDDIR)/mod_data.o $(BUILDDIR)/gto_eval.o $(BUILDDIR)/grid_gen.o
# Descendant submodules of exchange_cosx_impl (K / K_LR / K_SR / force
# builders, split out of the 4127-line parent). Each needs the PARENT's
# .smod, so each depends on the parent object, not just mod_integrals.o.
$(BUILDDIR)/mod_integrals_exchange_cosx_k.o     : $(BUILDDIR)/mod_integrals_exchange_cosx.o
$(BUILDDIR)/mod_integrals_exchange_cosx_lr.o    : $(BUILDDIR)/mod_integrals_exchange_cosx.o
$(BUILDDIR)/mod_integrals_exchange_cosx_sr.o    : $(BUILDDIR)/mod_integrals_exchange_cosx.o
$(BUILDDIR)/mod_integrals_exchange_cosx_force.o : $(BUILDDIR)/mod_integrals_exchange_cosx.o
$(BUILDDIR)/mod_integrals_force.o    : $(BUILDDIR)/mod_integrals.o $(BUILDDIR)/mod_data.o $(BUILDDIR)/mod_exchange.o
$(BUILDDIR)/mod_integrals_hessian.o  : $(BUILDDIR)/mod_integrals.o $(BUILDDIR)/mod_data.o
$(BUILDDIR)/hessian_1e.o             : $(BUILDDIR)/mod_integrals.o
$(BUILDDIR)/cphf.o                   : $(BUILDDIR)/mod_integrals.o $(BUILDDIR)/mod_data.o $(BUILDDIR)/mod_exchange.o
$(BUILDDIR)/mod_integrals_df_setup.o        : $(BUILDDIR)/mod_integrals.o $(BUILDDIR)/mod_mem_predict.o $(BUILDDIR)/mod_meminfo.o $(BUILDDIR)/mod_exchange.o $(BUILDDIR)/mod_data.o $(BUILDDIR)/mod_profile.o $(BUILDDIR)/mod_basis_files.o
$(BUILDDIR)/mod_integrals_df_coulomb.o      : $(BUILDDIR)/mod_integrals.o $(BUILDDIR)/mod_data.o $(BUILDDIR)/mod_profile.o
$(BUILDDIR)/mod_integrals_df_exchange.o     : $(BUILDDIR)/mod_integrals.o $(BUILDDIR)/mod_localization.o $(BUILDDIR)/mod_data.o $(BUILDDIR)/mod_profile.o
$(BUILDDIR)/mod_integrals_df_force_store.o  : $(BUILDDIR)/mod_integrals.o $(BUILDDIR)/mod_mem_predict.o $(BUILDDIR)/mod_meminfo.o $(BUILDDIR)/mod_data.o $(BUILDDIR)/mod_profile.o $(BUILDDIR)/mod_exchange.o
$(BUILDDIR)/mod_integrals_df_force_direct_j.o  : $(BUILDDIR)/mod_integrals.o $(BUILDDIR)/mod_data.o $(BUILDDIR)/mod_profile.o
$(BUILDDIR)/mod_integrals_df_force_direct_jk.o : $(BUILDDIR)/mod_integrals.o $(BUILDDIR)/mod_data.o $(BUILDDIR)/mod_profile.o
$(BUILDDIR)/mod_integrals_cosmo.o : $(BUILDDIR)/mod_integrals.o $(BUILDDIR)/mod_data.o
# mod_cosmo_cavity split (src/solvent/cavity/) - dependency chain:
# constants <- dual/radii/geodesic <- clip/gaubon <- cavity (orchestrator)
$(BUILDDIR)/mod_cosmo_dual.o     : $(BUILDDIR)/mod_cosmo_constants.o
$(BUILDDIR)/mod_cosmo_radii.o    : $(BUILDDIR)/mod_cosmo_constants.o
$(BUILDDIR)/mod_cosmo_geodesic.o : $(BUILDDIR)/mod_cosmo_constants.o
$(BUILDDIR)/mod_cosmo_clip.o     : $(BUILDDIR)/mod_cosmo_constants.o $(BUILDDIR)/mod_cosmo_dual.o
$(BUILDDIR)/mod_cosmo_gaubon.o   : $(BUILDDIR)/mod_cosmo_constants.o $(BUILDDIR)/mod_cosmo_dual.o
$(BUILDDIR)/mod_cosmo_cavity.o   : $(BUILDDIR)/mod_cosmo_constants.o $(BUILDDIR)/mod_cosmo_dual.o \
                                    $(BUILDDIR)/mod_cosmo_radii.o $(BUILDDIR)/mod_cosmo_geodesic.o \
                                    $(BUILDDIR)/mod_cosmo_clip.o $(BUILDDIR)/mod_cosmo_gaubon.o
# mod_smd_cds split (src/solvent/smd/) - dependency chain:
# tables <- sts/daareal <- cds (orchestrator)
$(BUILDDIR)/mod_smd_tables.o  : $(BUILDDIR)/mod_cosmo_cavity.o
$(BUILDDIR)/mod_smd_sts.o     : $(BUILDDIR)/mod_smd_tables.o
$(BUILDDIR)/mod_smd_daareal.o : $(BUILDDIR)/mod_smd_tables.o
$(BUILDDIR)/mod_smd_cds.o     : $(BUILDDIR)/mod_smd_tables.o $(BUILDDIR)/mod_smd_sts.o $(BUILDDIR)/mod_smd_daareal.o \
                                 $(BUILDDIR)/mod_cosmo_cavity.o
# mod_cosmo split (src/solvent/cosmo_impl/) - dependency chain:
# state <- init/scf/force <- cosmo (orchestrator); solvents is stateless
$(BUILDDIR)/mod_cosmo_solvents.o : $(BUILDDIR)/mod_data.o
$(BUILDDIR)/mod_cosmo_init.o  : $(BUILDDIR)/mod_integrals.o $(BUILDDIR)/mod_cosmo_cavity.o $(BUILDDIR)/mod_cosmo_state.o \
                                 $(BUILDDIR)/mod_data.o $(BUILDDIR)/mod_meminfo.o
$(BUILDDIR)/mod_cosmo_scf.o   : $(BUILDDIR)/mod_integrals.o $(BUILDDIR)/mod_cosmo_cavity.o $(BUILDDIR)/mod_cosmo_state.o
$(BUILDDIR)/mod_cosmo_force.o : $(BUILDDIR)/mod_integrals.o $(BUILDDIR)/mod_cosmo_cavity.o $(BUILDDIR)/mod_cosmo_state.o
$(BUILDDIR)/cosmo.o : $(BUILDDIR)/mod_cosmo_state.o $(BUILDDIR)/mod_cosmo_init.o $(BUILDDIR)/mod_cosmo_scf.o \
                       $(BUILDDIR)/mod_cosmo_force.o $(BUILDDIR)/mod_cosmo_solvents.o
$(BUILDDIR)/mod_exchange.o  : $(BUILDDIR)/mod_integrals.o
$(BUILDDIR)/mod_density_fitting.o : $(BUILDDIR)/mod_integrals.o
# mod_localization has no integrals/DF dependency at all (see its own
# header) - only needs mod_data.o for MOL_info's atoms(:)/nAtoms.
$(BUILDDIR)/mod_localization.o : $(BUILDDIR)/mod_data.o
$(BUILDDIR)/tool_io.o  : $(BUILDDIR)/mod_data.o
$(BUILDDIR)/mod_checkpoint.o : $(BUILDDIR)/mod_data.o
$(BUILDDIR)/mod_molden.o : $(BUILDDIR)/mod_data.o $(BUILDDIR)/mod_integrals.o
$(BUILDDIR)/guess.o    : $(BUILDDIR)/mod_data.o $(BUILDDIR)/mod_integrals.o
$(BUILDDIR)/SCF.o      : $(BUILDDIR)/mod_data.o $(BUILDDIR)/mod_integrals.o $(BUILDDIR)/mod_exchange.o $(BUILDDIR)/mod_density_fitting.o $(BUILDDIR)/mod_profile.o $(BUILDDIR)/mod_scf_history.o $(BUILDDIR)/mod_xc.o $(BUILDDIR)/mod_vv10.o
$(BUILDDIR)/scf_fock.o : $(BUILDDIR)/mod_data.o $(BUILDDIR)/mod_integrals.o $(BUILDDIR)/mod_exchange.o $(BUILDDIR)/mod_density_fitting.o $(BUILDDIR)/mod_profile.o $(BUILDDIR)/cosmo.o $(BUILDDIR)/mod_vv10.o
$(BUILDDIR)/solveKS.o  : $(BUILDDIR)/mod_data.o
$(BUILDDIR)/grid_gen.o : $(BUILDDIR)/mod_data.o
$(BUILDDIR)/grid_dynamic.o : $(BUILDDIR)/mod_data.o $(BUILDDIR)/mod_profile.o
$(BUILDDIR)/gto_eval.o : $(BUILDDIR)/mod_data.o $(BUILDDIR)/mod_meminfo.o $(BUILDDIR)/mod_integrals.o
$(BUILDDIR)/DFT.o      : $(BUILDDIR)/mod_data.o $(BUILDDIR)/mod_xc.o $(BUILDDIR)/mod_vv10.o $(BUILDDIR)/mod_vv10_aca.o $(BUILDDIR)/mod_profile.o
$(BUILDDIR)/DFT_vv10_pass.o : $(BUILDDIR)/mod_data.o $(BUILDDIR)/mod_vv10.o $(BUILDDIR)/mod_vv10_aca.o $(BUILDDIR)/mod_profile.o
$(BUILDDIR)/DFT_shellpair.o: $(BUILDDIR)/mod_data.o $(BUILDDIR)/mod_xc.o
$(BUILDDIR)/mod_xc.o  : $(BUILDDIR)/mod_vv10.o $(LIBXCDIR)/build/libxcf03.a
$(BUILDDIR)/force.o    : $(BUILDDIR)/mod_data.o $(BUILDDIR)/mod_integrals.o $(BUILDDIR)/mod_exchange.o $(BUILDDIR)/mod_profile.o $(BUILDDIR)/mod_vv10.o $(BUILDDIR)/mod_xc.o $(BUILDDIR)/cosmo.o
$(BUILDDIR)/engine.o   : $(BUILDDIR)/mod_data.o $(BUILDDIR)/mod_integrals.o $(BUILDDIR)/mod_xc.o $(BUILDDIR)/mod_exchange.o $(BUILDDIR)/mod_profile.o $(BUILDDIR)/mod_meminfo.o $(BUILDDIR)/mod_mem_predict.o $(BUILDDIR)/mod_checkpoint.o $(BUILDDIR)/mod_molden.o $(BUILDDIR)/mod_vv10.o $(BUILDDIR)/cosmo.o $(BUILDDIR)/cphf.o $(BUILDDIR)/mod_dispersion.o $(BUILDDIR)/mod_dispersion_d3.o $(BUILDDIR)/mod_dispersion_d4.o $(BUILDDIR)/mod_gcp.o $(BUILDDIR)/mod_smd_cds.o $(BUILDDIR)/mod_basis_files.o

# engine.f90 now contains only the library entry points (EngineUp,
# reset_engine_state) - no "program" statement - so the library objects
# get archived into libengine.a, and the example driver below supplies
# the program entry point.
$(BUILDDIR)/libengine.a: $(objects)
	ar rcs $@ $(objects)

# mod_engine_input* (src/io): input-file parsing, split out of
# run_engine.f90 so a future multi-step driver (OPT/BOMD) can reuse the
# parser without re-implementing it - see mod_engine_input.f90's header.
$(BUILDDIR)/mod_engine_input_atoms_namelist.o: mod_engine_input_atoms_namelist.f90 $(BUILDDIR)/mod_engine_input_elements.o | $(BUILDDIR)
	$(FORT90) $(tagCOMP) $< -o $@

$(BUILDDIR)/mod_engine_input_atoms_xyz.o: mod_engine_input_atoms_xyz.f90 $(BUILDDIR)/mod_engine_input_elements.o | $(BUILDDIR)
	$(FORT90) $(tagCOMP) $< -o $@

$(BUILDDIR)/mod_engine_input_atoms_column.o: mod_engine_input_atoms_column.f90 $(BUILDDIR)/mod_engine_input_elements.o | $(BUILDDIR)
	$(FORT90) $(tagCOMP) $< -o $@

$(BUILDDIR)/mod_engine_input.o: mod_engine_input.f90 $(BUILDDIR)/mod_engine_input_types.o $(BUILDDIR)/mod_engine_input_block.o \
                                 $(BUILDDIR)/mod_engine_input_elements.o $(BUILDDIR)/mod_engine_input_atoms_namelist.o \
                                 $(BUILDDIR)/mod_engine_input_atoms_xyz.o $(BUILDDIR)/mod_engine_input_atoms_column.o \
                                 $(BUILDDIR)/cosmo.o | $(BUILDDIR)
	$(FORT90) $(tagCOMP) $< -o $@

$(BUILDDIR)/run_engine.o: run_engine.f90 $(BUILDDIR)/mod_integrals.o $(BUILDDIR)/mod_profile.o $(BUILDDIR)/mod_scf_history.o $(BUILDDIR)/mod_version.o $(BUILDDIR)/mod_vv10.o $(BUILDDIR)/cosmo.o \
                          $(BUILDDIR)/mod_engine_input.o $(BUILDDIR)/mod_engine_input_types.o $(BUILDDIR)/mod_engine_input_block.o | $(BUILDDIR)
	$(FORT90) $(tagCOMP) $< -o $@


Engine : $(BUILDDIR)/run_engine.o $(io_objects) $(BUILDDIR)/mod_version.o $(BUILDDIR)/libengine.a $(LIBDIR)/libcint_ecp.a $(LIBCINTDIR)/build/libcint.a $(OPENBLASDIR)/libopenblas.a $(LIBXCDIR)/build/libxcf03.a $(TOMLF_PREREQ)
	$(FORT90) -o Engine -fbacktrace -fopenmp $(BUILDDIR)/run_engine.o $(io_objects) $(BUILDDIR)/mod_version.o $(BUILDDIR)/libengine.a $(LIBCINTDIR)/build/libcint.a $(LIBDIR)/libcint_ecp.a $(OPENBLASDIR)/libopenblas.a $(LIBXCDIR)/build/libxcf03.a $(LIBXCDIR)/build/libxc.a -lpthread -L$(DFTD3DIR)/build -Wl,-rpath,'$$ORIGIN/$(DFTD3DIR)/build' -ls-dftd3 -L$(MCTCDIR)/build -lmctc-lib -L$(GCPDIR)/build -Wl,-rpath,'$$ORIGIN/$(GCPDIR)/build' -lgcp -L$(DFTD4DIR)/build -Wl,-rpath,'$$ORIGIN/$(DFTD4DIR)/build' -ldftd4 -L$(THIRDPARTY_INSTALL)/lib -lmulticharge $(TOMLF_LINKFLAGS)

.PHONY : clean
clean :
	rm -f Engine $(objects) $(io_objects) $(BUILDDIR)/run_engine.o $(BUILDDIR)/libengine.a $(BUILDDIR)/*.mod \
	      $(BUILDDIR)/mod_version.f90 $(BUILDDIR)/mod_version.o \
	      $(BUILDDIR)/nr_ecp.o $(BUILDDIR)/nr_ecp_deriv.o $(BUILDDIR)/ecp2f.o $(LIBDIR)/libcint_ecp.a
