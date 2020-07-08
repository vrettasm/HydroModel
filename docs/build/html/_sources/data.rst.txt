Data Format
===========
The input data are given in a 'csv' file format, with the following structure:

-----

.. csv-table::
   :header: " ID ", " Date ", " Precipitation ", " Water_table_depths "
   :widths: auto
   :align: center

   1, 733682, 0, -4.6583
   2, 733683, 0.2, -4.6583
   3, 733684, 0.35, -4.6626
   ..., ..., ..., ...

-----

The **Date** column is a 'datenum' object. The precipitation is reported in [L:cm]
and the water table depths in [L:m], increasing downwards. The negative sign means
its underground, but it is removed in the code so it doesn't play a role here. The
datafile can be given in two ways:

    1. through the command line using the '--data <datafile>' option, or
    2. through the 'input_parameters.json' file using:
        * "Data_Filename": "path/to/the/datafile.csv"

If by accident the user passes two files (one with either option); the process will
use the datafile given through the command line (option 1).

.. note::
   The *precipitation* values are not allowed to contain NaN. If this is the case
   the process will terminate with an Error.

.. warning::
   The *datenum* object here is inherited from the MATLAB datafiles, which means
   that they had to be converted in the Python code, using pandas.to_datetime(),
   to match the same dates. This is happening in the *Simulation.setupModel()*
   method with: ::

    timestamps = pd.to_datetime(r_datenum - 719529, unit='D')

   The value *719529* is MATLAB's 'datenum' value of the "Unix epoch" starting
   (1970-01-01), which is the default origin for pandas.to_datetime().