# HydroModel - Berkeley
---
Hydrological Model (version 1.0.0)

This project implements the Python code of the underground (stochastic)
hydrological model that was developed during my postdoctoral tenure at
the Dept. of Earth & Planetary Science, U. C. Berkeley, (2013 - 2016).

There might be updates in the future, but this first version is now fully
operational. The data (water table depths, precipitation values, etc.) are
not available online and could be accessible only through communication
to Prof. Inez Fung (PI of the project).

## Data Format
---
The input data must be provided in a csv-file format with the following structure:

|   ID  |  Date  |  Precipitation  |  Water Table Depth  |
| :---: | :----: | :-------------: | :-----------------: |
| 1     | 733681 |          0.00   |            -10.4001 |
| 2     | 733682 |          0.01   |            -10.4151 |
| 3     | 733683 |          0.24   |            -10.4151 |
| ...   |  ...   |    ...          |               ...   |

The *Date* is a 'datenum' object. The *Precipitation* is given in [L: cm] and the
*Water Table Depths* are in [L: m] units (increasing downwards). The negative sign
indicates underground values but in the code is removed.

**Note:**
    We need to convert the dates (from MATLAB to Python). The value `719529` is
    MATLAB's datenum value of the "Unix epoch" start (1970-01-01), which is the
    default origin for pandas.to_datetime(). Hence:

    timestamps = pd.to_datetime(r_datenum - 719529, unit='D')

**Warning:**
   The precipitation column is not allowed to have NaN values. In such case an
   error will be raised and terminate the program.

## Installation
---
There are two options to install the software.

1. The easiest way is to visit the GitHub web-page of the project and
[download the code](https://github.com/vrettasm/HydroModel/archive/master.zip)
in zip format. This option does not require a prior installation of git
on the computer.

2. Alternatively one can clone the project directly using git as follows:

    `$ git clone https://github.com/vrettasm/HydroModel.git`

### Required packages

The recommended version is Python version 3.6+.
Other required packages are:

>
> Numpy, Scipy, Numba, Pandas, h5py, json
>

## How to run
---
To execute the program, first navigate to the main directory of the project
(i.e. where the berkeley_hydro_main.py is located), and then run the following
command:

    $ python3 -O berkeley_hydro_main.py --params ./model_parameters/input_parameters.json

The ‘-O’ option is not necessary, but it could speed up the execution of the simulation.

**Note:**
This assumes that the ‘input_parameters.json’ file includes the datafile location in the:

    “Data_Filename”: “path/to/datafile.csv”

## References
---
The work is described with details in two (open access) publications:

1. Michail D. Vrettas and Inez Y. Fung (2015). "Toward a new parameterization
of hydraulic conductivity in climate models: Simulation of rapid groundwater
fluctuations in Northern California".
Journal of Advances in Modeling Earth Systems (JAMES), vol. (7), issue (4),
pp:2105-2135, https://doi.org/10.1002/2015MS000516.

2. Michail D. Vrettas and Inez Y. Fung (2017). "Sensitivity of transpiration
to subsurface properties: Exploration with a 1-D model".
Journal of Advances in Modeling Earth Systems (JAMES), vol. (9), issue (2),
pp:1030-1045. https://doi.org/10.1002/2016MS000901.

### Contact
---
If you have any questions / comments please contact me at: vrettasm@gmail.com

M. Vrettas, PhD.


```
{
  "firstName": "John",
  "lastName": "Smith",
  "age": 25
}
``` 