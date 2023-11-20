#!/usr/bin/env -S bash -l

export MAMBA_HOME="/scratch/fl3/ycoles/fp2-wheat/containers/mamba"
export MAMBA_ROOT_PREFIX="${MAMBA_HOME}/prefixes/pawsey"

export PATH="${MAMBA_HOME}/prefixes/pawsey/bin:$PATH"
source "${MAMBA_HOME}/prefixes/pawsey/etc/profile.d/micromamba.sh"

micromamba activate "${MAMBA_HOME}/prefixes/pipeline"
