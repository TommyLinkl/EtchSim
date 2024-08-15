# Semiconductor Etching Simulation

## Overview

This project focuses on simulating the degradation and etching processes of wurtzite (wz)-InP nanoplatelets (NPLs) in solution using kinetic Monte Carlo simulations. The work involves understanding both controlled (layer-by-layer) and uncontrolled etching processes and their impact on the optoelectronic properties of the NPLs.


## Code Usage

### Prerequisites

To run the simulation, you need Python 3.x and some libraries. You can install the necessary libraries using the provided `requirements.txt` file.

```
pip install -r requirements.txt
```

### Running the Simulation

1. **Prepare Input File**: Ensure `input.json` is in the specified directory.
2. **Execute**: Run the script using:
   ```
   python main.py CALC_DIR/
   ```
3. **Post-processing**: If `process_stats_now` is set to 0, then post-processing is needed to generate the statistics. Run the following with the same `CALC_DIR/`
   ```
   python postproc.py CALC_DIR/
   ```

### Input File Description
The simulation requires an input JSON file (`input.json`) that specifies various settings and parameters. Below is a description of the expected structure and key fields of the input file: 

- **`calc_setting`**: Controls the simulation's logging and execution settings.
  - **`verbosity`**: Level of output detail (0 for minimal, higher values for more details).
  - **`runtime_flag`**: Indicates whether to output runtime. 
  - **`write_every`**: Frequency (in steps) to write intermediate results to the trajectories. 
  - **`process_stats_now`**: Flag to determine if statistics should be processed immediately (1 to process, 0 to delay and do postprocessing).
  - **`random_seed`**: Seed for random number generation in the kinetic Monte Carlo process to ensure reproducibility.

- **`NPL_setting`**: Specifies the path to the file containing lattice site data.
  - **`read_sites_from`**: Path to the `.pkl` file that contains the site information (e.g., `"../sites.pkl"`).

- **`sim_params`**: Contains parameters for the simulation. Please see the section below for parameter details. 
  - **`epsilon`** ($\epsilon$): Bonding energy parameter, usually set to 1 as the unit of energy. Estimated range: 0.1-0.5 eV.
  - **`mu_In`**, **`mu_P`** ($\mu_{In}$, $\mu_{P}$): Chemical potentials for Indium and Phosphorus atoms.Tested between -4 and -1.
  - **`T`**: Temperature in energy units. Set between 0.05 and 0.25 eV.
  - **`max_steps`**: Maximum number of steps for the simulation.


## Scientific Goals

Experimentally, the NPLs were prepared in a liquid cell. Via high-resolution TEM measurements, two etching processes are observed: first is a fast, uncontrolled etching at high electron dose rate, and second, a more common layer-by-layer etching trend in InP nanoplatelets. 

Theoretical simulation is used to support the observed etching mechanism. Further, it has been well established that surface has a big effect on nanosystems' optoelectronics properties. Understanding layer-by-layer etching vs. rough (irregular/uncontrolled) etching in wz-InP NPLs provides an additional control over their optoelectronic properties. The nanoplatelet size used in these simulations are: 15.4 nm in diameter, 3.5 nm in thickness. 

## Simulation Details

The simulation uses the Gillespie Algorithm. At each step, the surface semiconductor atoms' detachment and attachment rates are calculated. The total rate and cumulative sum are calculated, which were then used to determine the time increment and select events randomly. The lattice information is updated after each iteration. This is then repeated until no atoms (or no more detachment or attachment events are available). 

### Performance

- **Runtime**: ~40 milliseconds per KMC iteration on AMD EPYC 7763 (Milan) CPUs.
- **Memory Utilization**: ~400 MB + 8 MB per thousand iterations.
- **Data Handling**: Storage-efficient compressed files are used. In addition, the user has the option to choose real-time processing or post-processing of the statistics.
