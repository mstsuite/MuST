#!/bin/bash

# Define the lattice constants range
for lattice in $(seq 5.50 0.10 7.00); do
    lattice_dir="${lattice}"  # Folder for this lattice constant
    echo "Setting up and submitting job for lattice constant: $lattice"

    # Create a directory for this lattice constant if it doesn't exist
    mkdir -p "$lattice_dir"

    # Copy necessary files into the new directory
    cp TiVNbTa_mt_v i_mst position.dat aengel_slurm_cpa "$lattice_dir"

    # Modify position.dat in the new directory
    sed -i "1s/.*/$lattice/" "$lattice_dir/position.dat"

    # Submit the job from the new directory
    cd "$lattice_dir"
    sbatch aengel_slurm_cpa
    cd ..

done

echo "All jobs submitted."

# Reminder to run generate_input_file.py after jobs complete
echo "Once all jobs are finished, run: python3 generate_input_file.py"

