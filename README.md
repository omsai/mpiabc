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

## Installation

Install `spack` and the MPI-ABC dependencies into a new spack environment
called "mpi-abc":\n

```bash
git clone -c feature.manyFiles=true https://github.com/spack/spack.git
cd spack
source share/spack/setup-env.sh
spack compiler find
spack env create mpi-abc
spack env activate mpi-abc
spack add mpich
spack add gsl+external-cblas
spack add dmtcp
spack add doxygen
spack add ninja
spack concretize
spack install
```

Then run `ninja` to build both the executable `model` and this PDF
documentation `latex/refman.pdf`:\n

```bash
ninja
```

## Usage

The included example Lotka-Volterra model illustrates the 4 steps necessary to
integrate a model into the MPI-ABC framework:

1. Fit priors matching parameters in the model to the calibrated.
2. Input samples into the model.
3. Run the model.
4. Convert model output into one-or-more distances; the distances are also
   called summary statistics.
   
Run the `model` generated by `ninja` with:\n

```
./model
```
