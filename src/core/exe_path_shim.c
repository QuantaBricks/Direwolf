/*===============================
 * Wrapper around readlink("/proc/self/exe", ...), called from
 * mod_basis_files.f90 to find data/bases/ relative to where the
 * running Engine binary actually lives, not the caller's current
 * working directory - argv[0] alone can't do this reliably (it may be
 * a bare "Engine" found via PATH, or a relative path like "../../Engine"
 * that's only meaningful from the shell's own cwd at invocation time),
 * but /proc/self/exe is always resolved by the kernel to an absolute,
 * symlink-free path regardless of how the process was launched.
 *
 * Returns the path length on success (buf is NUL-terminated), or -1 on
 * failure (e.g. non-Linux /proc, though this codebase already targets
 * Linux elsewhere) - the Fortran caller falls back to the historical
 * CWD-relative "data/bases/" in that case.
 *===============================*/
#include <unistd.h>

long engine_get_exe_path(char *buf, long bufsize) {
    ssize_t n = readlink("/proc/self/exe", buf, (size_t)(bufsize - 1));
    if (n < 0) return -1;
    buf[n] = '\0';
    return (long)n;
}
