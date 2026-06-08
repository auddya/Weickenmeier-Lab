
# Periventricular White Matter Hyperintensity Simulations Base Code

This repository contains the base configuration and source files required to run periventricular white matter hyperintensity (PVWMH) simulations using **Abaqus 2022**. The computational model integrates ventricular growth kinematics with hyperelastic constitutive modeling.

## Repository Contents

| File | Description |
| --- | --- |
| `BrainFull.inp` | Main segmentation and mesh file (`.inp`) generated using Simpleware. |
| `UMAT_new.for` | Fortran user material (UMAT) subroutine containing separate functions for implementing ventricular growth and the hyperelastic constitutive models. |
| `initial-area.inp` | Precomputed initial area definitions for the thin shells representing the ventricles. |
| `script.sh` | Slurm batch script for submitting the simulation to a high-performance computing (HPC) cluster. |

## Prerequisites

To run these simulations, the following software and environments are required:

* **Abaqus 2022** (Standard)
* A compatible **Fortran Compiler** (e.g., Intel Fortran) linked to Abaqus for compiling the user subroutine.
* **Slurm Workload Manager** (if executing via the provided cluster script).

## Usage

### Submitting to an HPC Cluster

To run the simulation on a cluster environment, submit the provided batch script using `sbatch`:

```bash
sbatch script.sh

```

*(Note: You may need to verify that the module loads and Abaqus environment variables inside `script.sh` are configured correctly for your specific cluster's architecture.)*

### Running Locally (Command Line)

If you are testing the UMAT or running a scaled-down version of the simulation locally, use the standard Abaqus command line execution:

```bash
abaqus job=BrainFull user=UMAT_new.for interactive

```

### How to use this:
You can simply copy the text inside the bounds above and save it as `README.md` in the root of your `master` branch. GitHub will automatically render the formatting, icons, and code blocks.
