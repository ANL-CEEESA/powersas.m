# Installation

### 1. System requirements
Matlab (7.1 or later) or GNU Octave (4.0.0 or later).

### 2. Installation
 * Extract source code to a directory.
 * Enter the directory in Matlab or GNU Octave.
 * Execute command `setup`. You will see the following sub-directories:
    * `/data`: Stores test system data, simulation settings data, etc.
    * `/example`: Some examples of using PowerSAS.m.
    * `/output`: Stores test result data.
    * `/internal`: Internal functions of PowerSAS.m computation core.
    * `/util`: Auxiliary functions including data I/O, plotting, data conversion, etc.
    * `/logging`: Built-in logging system.
    * `/doc`: Documentation.

### 3. Test
 * Execute command `initpowersas` to initialize the environment, then execute `test_powersas` to run tests. You should expect all tests to pass.

### 4. Initialization
To initialize PowerSAS.m, add the main directory of PowerSAS.m to your Matlab or GNU Octave path and run the command `initpowersas`. This will ensure that all the functions of PowerSAS.m are added to the path and thus callable.
