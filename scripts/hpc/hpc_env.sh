# Sourced by the HPC scripts. Names the interpreter and compiler; the default batch python3 is never trusted.
# It only reads the shared environment: nothing here activates, installs into or changes it.
#
#   HPC_TOOLCHAIN       environment providing Python and GCC (default: /apps/anaconda3/envs/gcc)
#   HPC_PYTHON          base interpreter for .hpc-venv (default: $HPC_TOOLCHAIN/bin/python)
#   HPC_PYTHON_VERSION  version it must report (default: 3.11, the wheels in wheelhouse/py311)
#   MODULES             extra modules, loaded later so they take priority, e.g. a GCC 12.3.0 module

HPC_TOOLCHAIN="${HPC_TOOLCHAIN:-/apps/anaconda3/envs/gcc}"
HPC_PYTHON="${HPC_PYTHON:-$HPC_TOOLCHAIN/bin/python}"
HPC_PYTHON_VERSION="${HPC_PYTHON_VERSION:-3.11}"
if [[ -d "$HPC_TOOLCHAIN/bin" ]]; then
  PATH="$HPC_TOOLCHAIN/bin:$PATH"
  export PATH
fi
