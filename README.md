# Overview

Approximate Bayesian computing (ABC) with sequential Monte Carlo (SMC) sampling
is an increasingly popular and principled approach for "tuning" slow-running
stochastic models with a large number of parameters to multiple datasets.  This
process of choosing intervals of parameter values to match the datasets is also
sometimes referred to as "model calibration".

Existing higher-level language implementations of ABC-SMC samplers lack a few
desirable features:

1. *Debugging crashes* of samplers written in higher-level languages that occur
   several days of walltime is more challenging than using rich coredumps from
   compiled languages like C and C++.
2. *Checkpointing and resuming* calibrations is generally not supported by
   higher-level language ABC-SMC samplers.
3. *Parallel logging* facilities are also not supported.
4. *Multivariate parameters* like Dirichlet distributions e.g. of observed
   infection outcomes.

Such limitations jusfify the additional time and complexity involved with using
lower-level languages to support these and other desirable features with more
intricate control over implementation.

## Configure dependencies

### Polaris cluster

Load the modules for the MPI-ABC dependencies:\n

```bash
module use /soft/modulefiles
module load cray-python spack-pe-base ninja math_libs/gsl doxygen
```
### Local

If not using Polaris, fallback to a `spack` environment called "mpi-abc":\n

```bash
git clone -c feature.manyFiles=true https://github.com/spack/spack.git
cd spack
source share/spack/setup-env.sh
spack compiler find
spack env create mpi-abc
spack env activate mpi-abc
spack add mpich
spack add gsl+external-cblas
spack add doxygen
spack add ninja
# For tests
spack add py-gcovr
spack concretize
spack install
```

Then to use this `spack` environment in the future:\n

```shell
source spack/share/spack/setup-env.sh
spack env activate mpi-abc
```

## Build

Run `meson` to build both the executable `infer`:\n

```bash
BUILDDIR=../build-mpiabc
meson setup -Db_coverage=true $BUILDDIR
meson compile -C $BUILDDIR
# Optional locations of headers for editing source files with Emacs flycheck:
cp -a $BUILDDIR/.dir-locals.el .
```

## Documentation

Generate the Doxygen PDF at `$BUILDDIR/latex/refman.pdf`:\n

```bash
BUILDDIR=../build-mpiabc
meson compile -C $BUILDDIR docs
```

## Test

```bash
BUILDDIR=../build-mpiabc
meson test --suite mpiabc -C $BUILDDIR
```

Check test coverage with:\n

```bash
BUILDDIR=../build-mpiabc
meson setup -Db_coverage=true --reconfigure $BUILDDIR
meson test --suite mpiabc -C $BUILDDIR
ninja coverage-text -C $BUILDDIR &&
	cat $BUILDDIR/meson-logs/coverage.txt
```

Check for memory leaks with:\n

```bash
BUILDDIR=../build-mpiabc
CK_FORK=no valgrind --leak-check=full $BUILDDIR/test_abc
```

## API

The included example Lotka-Volterra model illustrates the 4 steps necessary to
integrate a model into the MPI-ABC framework:

1. Fit priors matching parameters in the model to the calibrated.
2. Fix the remaining parameters.
3. Check that the model produces outcome variables.
4. Convert model outcomes into one-or-more distances; the distances are also
   called summary statistics.
