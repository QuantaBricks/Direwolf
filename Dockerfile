# Works unmodified against both the full dev tree (this repo) and an
# exported release tree (release/Direwolf-<ver>/, see scripts/make_release.sh)
# - both share the same Makefile/third_party layout, so nothing here is
# tree-specific.
#
# Build: docker build -t direwolf .
# Run:   docker run --rm -v "$PWD":/data direwolf /data/input.inp /data/output.out

# ---- stage 1: build everything from vendored source -----------------------
FROM debian:bookworm-slim AS builder

RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential gfortran cmake git ca-certificates \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /src

# x86-64-v2 (SSE4.2/POPCNT, ~2009+), not the default -march=native: this
# binary has to run on whatever host later docker-runs the image, which is
# almost never the exact build machine (same reasoning as make_release.sh's
# binary-bundle build). Override per target arch, e.g.
# --build-arg MARCH=armv7-a+neon-vfpv4 for a 32-bit ARM (Raspberry Pi 2/3)
# image.
ARG MARCH=x86-64-v2

# The vendored third_party/ deps (OpenBLAS, libxc, libcint, toml-f,
# mctc-lib, dftd3/4, gcp, multicharge) are a self-contained ~30-min build
# (much longer under QEMU for a cross-arch image) that only changes when
# third_party/ itself does. Build it as its own cached layer, before the
# full source COPY, so iterating on src/ doesn't re-pay it every time.
COPY Makefile ./
COPY third_party ./third_party
RUN make MARCH="$MARCH" -j"$(nproc)" \
        third_party/libcint/build/libcint.a \
        third_party/OpenBLAS/libopenblas.a \
        third_party/libxc/build/libxcf03.a \
        third_party/toml-f/build/libtoml-f.a \
        third_party/mctc-lib/build/libmctc-lib.a \
        third_party/simple-dftd3/build/libs-dftd3.so \
        third_party/gcp/build/libgcp.so \
        third_party/multicharge/build/libmulticharge.a \
        third_party/dftd4/build/libdftd4.so

COPY . .
RUN make MARCH="$MARCH" -j"$(nproc)" Direwolf

# ---- stage 2: minimal runtime ----------------------------------------------
FROM debian:bookworm-slim AS runtime

# Direwolf links libgfortran/libgomp/libgcc_s dynamically (see Makefile's
# Direwolf: link line); everything else (libcint, OpenBLAS, libxc) is
# statically linked in, so this is the whole runtime dependency list - no
# gfortran/cmake/build-essential in the final image. libquadmath0 is not
# listed explicitly: on amd64 it comes in as a libgfortran5 dependency,
# and on armhf it does not exist at all (GCC ships no quad-math lib for
# 32-bit ARM) and is not needed - Direwolf uses no real(16).
RUN apt-get update && apt-get install -y --no-install-recommends \
        libgfortran5 libgomp1 libgcc-s1 \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

COPY --from=builder /src/Direwolf ./
COPY --from=builder /src/data ./data
COPY --from=builder /src/LICENSE ./
# The 3 LGPL runtime .so's (simple-dftd3, gcp, dftd4) - Direwolf's rpath is
# $ORIGIN-relative to these exact paths (see Makefile's -Wl,-rpath lines),
# so the directory layout below must be preserved as-is.
COPY --from=builder /src/third_party/simple-dftd3/build/libs-dftd3.so* ./third_party/simple-dftd3/build/
COPY --from=builder /src/third_party/gcp/build/libgcp.so* ./third_party/gcp/build/
COPY --from=builder /src/third_party/dftd4/build/libdftd4.so* ./third_party/dftd4/build/

ENTRYPOINT ["./Direwolf"]
