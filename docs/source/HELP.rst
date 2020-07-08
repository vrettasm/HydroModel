Help
====
This page provides some basic help to install and run the software.
The instructions are valid for Linux and MacOs operating systems.

Installation
------------
-   The easiest way to get the software is to visit the GitHub web-page
    of the project: https://github.com/vrettasm/HydroModel and download
    the code in zip format. Uncompress the file in any directory and you
    will have all the code in one place. This option does not require to
    have git installed on the computer.

-   Alternatively one can directly clone the project with git using: ::

    $ git clone https://github.com/vrettasm/HydroModel.git

Run the code
------------
From the main directory of the project (i.e. where the 'berkeley_hydro_main.py'
is located) run the following command, in the console: ::

    $ python3 -O berkeley_hydro_main.py --params ./model_parameters/input_parameters.json

The '-O' option is not necessary, but it could help to speed up the execution
of the simulation.

.. note::
   The 'input_parameters.json' file should contain all the model parameters for the simulation
   including the location (path) of the datafile.

Contact
-------
If you have any questions please contact me at: vrettasm@gmail.com