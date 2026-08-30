# Works unmodified against both the full dev tree (this repo) and an
# exported release tree (release/Engine-<ver>/, see scripts/make_release.sh)
# - both share the same Makefile/third_party layout, so nothing here is
# tree-specific.
#
# Build: docker build -t engine .
# Run:   docker run --rm -v "$PWD":/data engine /data/input.inp /data/output.out

# ---- stage 1: build everything from vendored source -----------------------
FROM debian:bookworm-slim AS builder

RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential gfortran cmake git ca-certificates \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /src
COPY . .

# x86-64-v2 (SSE4.2/POPCNT, ~2009+), not the default -march=native: this
# binary has to run on whatever host later docker-runs the image, which is
# almost never the exact build machine (same reasoning as make_release.sh's
# binary-bundle build).
RUN make MARCH=x86-64-v2 -j"$(nproc)" Engine

# ---- stage 2: minimal runtime ----------------------------------------------
FROM debian:bookworm-slim AS runtime

# Engine links libgfortran/libgomp/libquadmath/libgcc_s dynamically (see
# Makefile's Engine: link line); everything else (libcint, OpenBLAS, libxc)
# is statically linked in, so this is the whole runtime dependency list -
# no gfortran/cmake/build-essential in the final image.
RUN apt-get update && apt-get install -y --no-install-recommends \
        libgfortran5 libgomp1 libquadmath0 libgcc-s1 \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

COPY --from=builder /src/Engine ./
COPY --from=builder /src/data ./data
COPY --from=builder /src/LICENSE ./
# The 3 LGPL runtime .so's (simple-dftd3, gcp, dftd4) - Engine's rpath is
# $ORIGIN-relative to these exact paths (see Makefile's -Wl,-rpath lines),
# so the directory layout below must be preserved as-is.
COPY --from=builder /src/third_party/simple-dftd3/build/libs-dftd3.so* ./third_party/simple-dftd3/build/
COPY --from=builder /src/third_party/gcp/build/libgcp.so* ./third_party/gcp/build/
COPY --from=builder /src/third_party/dftd4/build/libdftd4.so* ./third_party/dftd4/build/

ENTRYPOINT ["./Engine"]
