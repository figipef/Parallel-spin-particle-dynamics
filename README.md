# Particle Dynamics Simulation

A high-performance simulation of a charged particle under electromagnetic fields.

## Table of Contents
- [Installation](#installation)
- [Usage](#usage)
- [Configuration](#configuration)  
  - [General Simulation Parameters](#general-simulation-parameters)  
  - [Particle Creation Parameters](#particle-creation-parameters)  
    - [Parameter Description](#parameter-description)  
    - [Parameter Options](#parameter-options)  
    - [Distribution Sizes](#distribution-sizes)  
  - [Diagnostics Parameters](#diagnostics-parameters)  
    - [Parameter Description](#parameter-description-1)  
    - [Parameter Options](#parameter-options-1)  
  - [Laser Parameters](#laser-parameters)  
    - [Parameter Description](#parameter-description-2)  
    - [Parameter Options](#parameter-options-2)  
- [Python Files](#python-files)
  - [File Overview](#file-overview)
  - [Dependencies](#dependecies)
- [Authors](#authors)

## Installation

1. Clone the repository:
   ```sh
   git clone https://github.com/figipef/Parallel-spin-particle-dynamics.git
   cd Parallel-spin-particle-dynamics
   cd build
   ```
2. Compile the project:
   ```sh
   make
   ```
3. Clean project folders:
   ```sh
   make clean
   ``` 

## Usage

Run the simulation using (in this case 8 processes):
```sh
./simulation_mpi.exe
``` 

Verbose options are `0`, `1`, `2`, `3` according to the user desire. It is `0` by default.
To apply verbosity run, where the number after the executable defines de verbosity option:
```sh
./simulation_mpi.exe 1
``` 

To run the python file:
```sh
cd ..
python main.py
```

## Configuration

Edit `build/input.txt` to modify the simulation parameters

### General Simulation Parameters

`TIME_STEP` is the normalized time between iterations
`TOTAL_TIME` is the normalized duration of the simulation
`STEPS_DIAG` is the number of timesteps between the writing of diagnostics (if it is 0 then no recording will be done)
`RADIATION_REACTION` is a `0` or a `1` to select if the user wants to take into account the radiation reaction effects
`FOLLOW_PARTICLE` is a `0` or a `1` to select if the user wants to get the position, momentum and spin of particles during the whole simulation
`RANDOM_PARTICLE` is a `0` or a `1` to select if the user wants to follow a predetermined one or random particles, respectively.
`PARTICLE_NUMBER` is either the particle index that the user wants to follow, or the number of random particles that should be followed (MAX 100)

### Particle Creation Parameters

#### Parameter Description

`NUMBER_OF_PARTICLES` is the total number of particles in the simulation
`*_DIST_TYPE` is the distribution type for a given `*` variable
`*_DIST_SIZE` is the size of the distribution in a single direction for a given `*` variable
`*_PREF_DIR` is the preferencial direction for static initalization of a given `*` variable

`*` is the identification of the variable in this case in can be:
`POSITION` | `MOMENTUM` | `SPIN`

#### Parameter Options

`POSITION_DIST_TYPE`:
 0 -> Uniform distribution
 1 -> Gaussian distribution
 2 -> Static (`POSITION_PREF_DIR`) 

`MOMENTUM_DIST_TYPE`:
 0 -> Uniform distribution
 1 -> Gaussian distribution
 2 -> Static (`MOMENTUM_PREF_DIR`)

`SPIN_DIST_TYPE`:(ALL WITH UNIT NORM)
 0 -> Uniform distribution 
 1 -> Static (`SPIN_PREF_DIR`) 
 2 -> Von Miser-Fisher (Central direction `SPIN_PREF_DIR`; Angular spread parameter, kappa `SPIN_DIST_SIZE`)

`*_PREF_DIR`:
 Should be written like: `X,Y,Z` where X,Y and Z can be any double.

#### Distribution Sizes

1. Uniform distribution:
	`*_DIST_SIZE` corresponds to the uniform distributuion size
2. Gaussian distribution:
	`*_DIST_SIZE` corresponds to the standard deviation. No particles are created with a `*` variable over 3 standard deviations.
3. Von Miser-Fisher:
	`*_DIST_SIZE` corresponds to the concentration parameter, kappa, of this distribution. Greater kappa means random spin direction are closer to the central direction.

### Diagnostics Parameters

Diagnostics are limited up to 2 to generate histograms in txt files. Values are not recorded if it exceeds the `BIN_MAX_` or the `BIN_MIN_`

#### Parameter Description

`PAR*` is the variable that want to be recorded
`BIN_NUMBER_*` is the number of bins for that variable
`BIN_MAX_*` is the largest value up to which the bins go to
`BIN_MIN_*` is the smallest value up to which the bins go to
`SHOW_E` and `SHOW_B` are booleans to either record the electric field and/or the magnetic field
`FIELD_BIN` is the size of the bins for the recorded fields
`FIELD_BIN_MAX` is the largest value up to which the field bins go to 
`FIELD_BIN_MIN` is the smallest value up to which the field bins go to 
`FIELD_BIN_DIR` is the direction on which the fields will be recorded

`*` is an integer that numbers the parameters. Must start at 1.

#### Parameter Options

`PAR*`:

 `p**` to record the position on the `**` direction
 `m**` to record the momentum on the `**` direction
 `s**` to record the spin on the `**` direction

 where `**` is in Cartesian Coordinates:
  1 -> x
  2 -> y
  3 -> z

`BIN_NUMBER_*`  is integer
`SHOW_E` and `SHOW_B` are booleans

The other parameters should be inputed as either integers or booleans

### Lasers Parameters

Lasers are limited to 10 (why would you ever need more);
All oscillatory fields are calculated with a cossine;
They are initialized at 0,0,0 ;
Particles will always be created in the direction of the longest laser in case of a finite pulse.

#### Parameter Description

`L*` is the type of laser/fields used
`K*` is the wave vector of the laser/fields used
`E*` is the inital electric field strength in normalized units
`B*` is the inital magnetic field strength in normalized units
`FREQ*` is the laser/fields frequency
`ENV_L*` is the length of the envelope of a finite laser pulse
`PHASE*` is multiple of PI/2 that is added as a phase

`*` is an integer that numbers the lasers. Must start at 1.

#### Parameter Options

`L*`:

 0 -> Constant magnetic and electric fields
 Initializion:

  `E*` + `K*` | `E*` + `B*`

 1 -> Oscillating magnetic and electric fields
 Initializion:
  
  `E*` + `K*` | `E*` + `B*` + `FREQ*` 

 2 -> Laser
 Initializion:
  
  `E*` + `K*`

 3 -> Laser in an envelope
 Initializion:
  
  `E*` + `K*` + `ENV_L*`

`E*` & `K*` & `B*`:
 Should be written like: `X,Y,Z` where X,Y,Z can be any double

`ENV_L*`:
 Can be either an integer or a double

`PHASE*`:
 Has to be an integer

## Python Files

The Python helper scripts are located in the `python_helper/` directory. These scripts assist with reading simulation data, performing analysis, and plotting results.

### File Overview

- `python_helper/Reader.py`  
  Responsible for reading simulation output files and input configuration files.

- `python_helper/Plotter.py`  
  Provides ready-to-use plotting functions for visualizing particle and spin dynamics.

- `main.py` (located in the root directory)  
  The main Python entry point that ties together reading, processing, and plotting — similar in spirit to a `main()` function in C++.

### Dependencies

If you encounter missing libraries, install them using:

```bash
pip install <library-name>
```

## Authors

This projected was made by:

Afonso santos
André Filipe
Rafael Fleming

Under the guidance of:

Profs Marija Vranic
Bernardo Barbosa
