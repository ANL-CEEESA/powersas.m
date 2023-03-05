# PowerSAS.m

### Documentation
[HTML](https://powersasm.readthedocs.io/en/latest/index.html)

[PDF](https://powersasm.readthedocs.io/_/downloads/en/latest/pdf/)

### What is PowerSAS.m?

**PowerSAS.m** is a robust, efficient and scalable power grid analysis framework based on semi-analytical solutions (SAS) technology. The **PowerSAS.m** is the version for MATLAB/Octave users. It currently provides the following functionalities (more coming soon!):

* **Steady-state analysis**, including power flow (PF), continuation power flow (CPF), contingency analysis.
* **Dynamic security analysis**, including voltage stability analysis, transient stability analysis, and flexible user-defined simulation.
* **Hybrid extended-term simulation** provides adaptive QSS-dynamic hybrid simulation in extended term with high accuracy and efficiency.

### Key features
* **High numerical robustness.** Backed by the SAS approach, the PowerSAS tool provides much better convergence than the tools using traditional Newton-type algebraic equation solvers when solving algebraic equations (AE)/ordinary differential equations (ODE)/differential-algebraic equations(DAE).  
* **Enhanced computational performance.** Due to the analytical nature, PowerSAS provides model-adaptive high-accuracy approximation, which brings significantly extended effective range and much larger steps for steady-state/dynamic analysis. PowerSAS has been used to solve large-scale system cases with 200,000+ buses.
* **Customizable and extensible.** PowerSAS supports flexible customization of grid analysis scenarios, including complex event sequences in extended simulation term.  

### Installation
#### 1. System requirements
Matlab (7.1 or later) or GNU Octave (4.0.0 or later).

#### 2. Installation
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

#### 3. Test
 * Execute command `initpowersas` to initialize the environment, then execute `test_powersas` to run tests. You should expect all tests to pass.

#### 4. Initialization
To initialize PowerSAS.m, add the main directory of PowerSAS.m to your Matlab or GNU Octave path and run the command `initpowersas`. This will ensure that all the functions of PowerSAS.m are added to the path and thus callable.