Data and Models
===============

1. Supported data formats
~~~~~~~~~~~~~~~~~~~~~~~~~

Currently PowerSAS.m supports extended PSAT (Matlab) data format.

In addition, PSS/E data format can be easily converted to the PSAT data format 
using the ``psse2mpc`` function by MATPOWER and the ``matpower2psat`` function
by PSAT. A sample code to convert a PSS/E *raw file to the PSAT format is below:
.. code:: matlab

    casename = 'ieee14'
    [mpc, warnings] = psse2mpc(horzcat(casename,'.raw'),horzcat(casename,'.m'));
    check = matpower2psat(horzcat(casename,'.m'), pwd);



Support for other formats and data format conversion features will be
added in future versions.

2. Extension of PSAT (Matlab) format
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

2.1 Automatic generation control (AGC) model
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Here PowerSAS provides a simple AGC model. It is named as ``Agc.con`` in
data files and it is a :math:`\text{N}\times \text{4}` matrix. Each
column is defined as below: Column \| Content ——-\| ————- 1 \| Bus index
2 \| Reciprocal of turbine governor gain on bus 3 \| Effective damping
ratio on bus 4 \| Reciprocal of AGC control time constant
