# MuST Silicon Full-Potential Benchmark

A MuST benchmark for crystalline silicon.

## Benchmark settings

| Setting | Value |
|---|---|
| System | Silicon (Si) |
| Method | Full potential (FP) |
| Local interaction zone (LIZ) | 338 atoms (default) |
| Magnetism | Spin-canted (noncollinear) |
| Energy contour | 32 energy points |
| SCF iterations | 5 |

## Contents

```
Si_nvagent/
├── i_mst                  # MuST main input file
├── ...                    # potential / position files used by i_mst
├── numerical_accuracy.sh  # extracts per-iteration SCF quantities from an o_ file
├── job.sh  # Job submission script (1 GH node 64 MPI task CPU run)
└── reference_run/         # reference result to compare against
    ├── o_n0000000_Si      # main output of the reference run
    └── k_n0000000_Si      # final SCF energy of the reference run
```

## Output files

| File | What it holds |
|---|---|
| `o_n0000000_Si` | Main output: every SCF iteration, including the Fermi energy, charges, moments and total energy |
| `k_n0000000_Si` | Final SCF energy |

## Running the benchmark

1. Run MuST with your usual job script, for example `ibrun mst2 < i_mst` inside a Slurm batch script.
2. When the run finishes, you should have `o_n0000000_Si` and `k_n0000000_Si` in the run directory.

## Checking numerical accuracy

`numerical_accuracy.sh` reads an `o_` file and prints these quantities for every SCF iteration:

- Fermi Energy DOS per Site
- Fermi Energy (Ryd)
- Integrated Electron Density in VP
- Integrated Moment Density in VP (x, y, z)
- Interstitial rho0
- Total Energy Per Atom

It also prints the final SCF energy from the matching `k_` file when that file is present.

```bash
chmod +x numerical_accuracy.sh

# Summary of your own run
./numerical_accuracy.sh o_n0000000_Si

# Summary of the reference run
./numerical_accuracy.sh reference_run/o_n0000000_Si

# Your run compared with the reference, iteration by iteration
./numerical_accuracy.sh o_n0000000_Si reference_run/o_n0000000_Si
```

Example output (first iteration of the reference run):

```
File: reference_run/o_n0000000_Si

=================== ITERATION 01 ===================
    Fermi Energy DOS per Site          =    2.780206969681
  Fermi Energy (Ryd)                   =     0.61761155968
    Integrated Electron Density in VP  =     3.80489480007
    Integrated Moment Density in VP    =     0.00000000000     0.00000000000    -0.00000000000
  Fermi Energy (Ryd)                   =     0.61761155968
        Total Energy Per Atom =    -576.75528820
...
```

The run goes on up to iteration 05.

In comparison mode, the script also prints a table of `|ΔE_tot|` and `|ΔE_F|` for each iteration against the reference.

## Notes

- With only 5 iterations the run is not expected to reach SCF convergence; the total energy still changes a lot from one iteration to the next. The benchmark measures **reproducibility**: compare each iteration with the same iteration of the reference run, not with a converged value.
- The full-potential case takes a long time on CPU, mainly in potential generation.
- Energies are in Rydberg.
