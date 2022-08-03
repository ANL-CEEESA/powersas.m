# Advanced Usage

### 1. Customize analysis with settings file/struct
PowerSAS.m lets you customize your simulation by providing a simulation settings interface to specify the events and scenarios to be analyzed. To use the customized simulation, call the `runPowerSAS` function as follows:
```Matlab
res=runPowerSAS('dyn',data,options,settings,snapshot)
``` 
Details are explained as follows:
#### 1.1 Settings file
The input argument `settings` can be a string specifying the settings file name, or a struct containing all the simulation settings.

When `settings` is a string, it should be a valid file name of an .m script file containing the settings. Some examples of the settings files can be found in the directory `/data`. The settings file should have the following variables:

* `eventList`: A gross list of events. ([more details](#variable-eventlist))
* `bsBus`: (for black start simulation only) Black start bus information. ([more details](#variable-bsbus))
* `bsSyn`: (for black start simulation only) Generator information on black start bus. ([more details](#variable-bssyn))
* `bsInd`: (for black start simulation only) Inductin motor on black start bus. ([more details](#variable-bsind))
* `Efstd`: Excitation potential of synchronous generators. ([more details](#variable-efstd))
* `evtLine`: List of line addition/outage events. ([more details](#variables-evtline-and-evtlinespec))
* `evtLineSpec`: Specifications of line addition/outage events. ([more details](#variables-evtline-and-evtlinespec))
* `evtZip`: List of static load addition/shedding events. ([more details](#variables-evtzip-evtzipspec-and-evtzipspec2))
* `evtZipSpec`: Specifications of static load addition/shedding events. ([more details](#variables-evtzip-evtzipspec-and-evtzipspec2))
* `evtZipSpec2`: Alternative specifications of static load addition/shedding events. ([more details](#variables-evtzip-evtzipspec-and-evtzipspec2))
* `evtInd`: List of induction motor addition/outage events. ([more details](#variables-evtind-and-evtindspec))
* `evtIndSpec`: Specifications of induction motor addition/outage events. ([more details](#variables-evtind-and-evtindspec))
* `evtSyn`: List of synchronous generator addition/outage events. ([more details](#variables-evtsyn-and-evtsynspec))
* `evtSynSpec`: Specifications of synchronous generator addition/outage events. ([more details](#variables-evtsyn-and-evtsynspec))
* `evtFault`: List of fault occurrence/clearing events. ([more details](#variables-evtfault-and-evtfaultspec))
* `evtFaultSpec`: Specifications of fault occurrence/clearning events. ([more details](#variables-evtfault-and-evtfaultspec))
* `evtDyn`: List of dynamic ramping events. ([more details](#variable-evtdyn))
* `evtDynPQ`: Specifications of PQ bus ramping. ([more details](#variable-evtdynpq))
* `evtDynPV`: Specifications of PV bus ramping. ([more details](#variable-evtdynpv))
* `evtDynInd`: Specifications of induction motor mechanical load ramping. ([more details](#variable-evtdynind))
* `evtDynZip`: Specifications of ZIP load ramping. ([more details](#variable-evtdynzip))
* `evtDynSh`: Specifications of shunt compensator ramping. ([more details](#variable-evtdynsh))
* `evtDynZipRamp`: Alternative specifications of ZIP load ramping. ([more details](#variable-evtdynramp))
* `evtDynTmech`: Specifications of generator mechanical torque ramping. ([more details](#variable-evtdyntmech))
* `evtDynPm`: Specifications of generator input active power ramping. ([more details](#variable-evtdynpm))
* `evtDynEf`: Specifications of generator excitation potential ramping. ([more details](#variable-evtdynef))
* `evtDynVref`: Specifications of exciter reference voltage ramping. ([more details](#variable-evtdynvref))
* `evtDynEq1`: Specifications of generator transient excitation potential ramping. ([more details](#variable-evtdyneq1))

#### 1.2 Settings struct
Alternatively, the `settings` can be a struct containing all the previous variables as its fields.


#### Example 1: Transient stability analysis (TSA) using settings file
In this example, we want to perform transient stability analysis on 2383-bus system. Here is the scenario of the TSA:
* The total simulation period is 10 seconds.
* At 0.5 s, apply faults on the starting terminals of lines 74 and 114, the fault resistance is 0.0 and the reactance is 0.02. At 0.75 s, clear the faults.
* At 1.5 s, apply fault on the starting terminal of line 1674, the fault resistance is 0.0 and the reactance is 0.1. At 1.95 s, clear the faults.

The settings file for this simulation is shown below. Other variables irrelevant to the fault events are omitted here for the sake of clarity.
```Matlab
eventList=[...
   1    0.0000  0.0000    0    1 0.0  0.0000
   1.1  0.5000  0.0000    6    1 0.0  0.0000
   1.2  0.7500  0.0000    6    2 0.0  0.0000
   1.3  1.5000  0.0000    6    3 0.0  0.0000
   1.4  1.9500  0.0000    6    4 0.0  0.0000
   3    10.00   0.0000   99    0 0.0  0.0000
];

% Fault event data
evtFault=[...
    1   1   2
    2   3   4
    3   5   5
    4   6   6
 ];

evtFaultSpec=[...
     114,  0.00, 0, 0.02,  0;
      74,  0.00, 0, 0.02,  0;
     114,  0.00, 0, 0.02,  1;
      74,  0.00, 0, 0.02,  1;
    1674,  0.00, 0, 0.1,   0;
    1674,  0.00, 0, 0.1,   1;
 ];
```

Assume the settings file is `settings_polilsh_tsa.m` and the system data file is `d_dcase2383wp_mod2_ind_zip_syn.m`. We can run the TSA as follows:
```Matlab
res_2383_st=runPowerSAS('pf','d_dcase2383wp_mod2_ind_zip_syn.m'); % Run steady-state
res_2383_tsa=runPowerSAS('dyn','d_dcase2383wp_mod2_ind_zip_syn.m',setOptions('hotStart',1),'settings_polilsh_tsa',res_2383_st.snapshot); % Hot start from existing steady-state

plotCurves(1,res_2383_tsa.t,res_2383_tsa.stateCurve,res_2383_tsa.SysDataBase,'v'); % plot the voltage magnitude curves
```

### 2. Extended-term simulation using hybrid QSS and dynamic engines
To accelerate computation — especially for extended-term simulation — PowerSAS.m provides an adaptive way to switch between QSS and dynamic engines in the course of a simulation. With this feature enabled, PowerSAS.m can switch to QSS simulation for better speed on detecting the fade-away of transients and switch back to dynamic simulation upon detecting transient events.

For more details on the technical approach, please refer to our paper:
* Hybrid QSS and Dynamic Extended-Term Simulation Based on Holomorphic Embedding, arXiv:2104.02877

Example 2 illustrates the use of PowerSAS.m to perform extended-term simulation.

#### Example 2: Extended-term simulation
We want to study the response of a 4-bus system under periodic disturbances. The entire simulated process is 500 seconds. Starting at 60 s and continuing until 270 s, the system undergoes events of adding/shedding loads every 30 s.  

The key settings of the simulation are:
```Matlab
% settings_d_004_2a_agc.m

eventList=[...
   1    0.0000  0.0000    0    1 0.0  0.0100
   6   60.0000  0.0000    2    1 0.0  0.0100
   7   90.0000  0.0000    2    2 0.0  0.0100
   8  120.0000  0.0000    2    3 0.0  0.0100
   9  150.0000  0.0000    2    4 0.0  0.0100
  10  180.0000  0.0000    2    1 0.0  0.0100
  11  210.0000  0.0000    2    2 0.0  0.0100
  12  240.0000  0.0000    2    3 0.0  0.0100
  13  270.0000  0.0000    2    4 0.0  0.0100
  18  500.0000  0.0000   99    0 0.0  0.0100
];

% Static load event data
evtZip=[...
   1    1    1    1
   2    1    2    2
   3    1    3    3
   4    1    4    4
];

evtZipSpec2=[...
   3 100.0000 100.0000 60.0000  0.0648  0.0648  0.0648  0.0359  0.0359  0.0359    0    1
   2 100.0000 100.0000 60.0000  0.0648  0.0648  0.0648  0.0359  0.0359  0.0359    0    1
   3 100.0000 100.0000 60.0000 -0.0648 -0.0648 -0.0648 -0.0359 -0.0359 -0.0359    0    1
   2 100.0000 100.0000 60.0000 -0.0648 -0.0648 -0.0648 -0.0359 -0.0359 -0.0359    0    1
];
```

First we run the simulation in full-dynamic mode and record time:
```Matlab
% Full dynamic simulation 
tagFullDynStart=tic;
res_004_fulldyn=runPowerSAS('dyn','d_004_2a_bs_agc.m',[]],'settings_d_004_2a_agc');
timeFullDyn=toc(tagFullDynStart);
```

Then we run the simulation in hybrid QSS & dynamic mode and record time:
```Matlab
% Hybrid simulation with dynamic-QSS switching
tagHybridStart=tic;
res_004=runPowerSAS('dyn','d_004_2a_bs_agc.m',setOptions('allowSteadyDynSwitch',1),'settings_d_004_2a_agc');
timeHybrid=toc(tagHybridStart);
```

Compare the results:
```Matlab
plotCurves(1,res_004_fulldyn.t,res_004_fulldyn.stateCurve,res_004_fulldyn.SysDataBase,'v');
plotCurves(2,res_004.t,res_004.stateCurve,res_004.SysDataBase,'v');
```

And compare the computation time:
```Matlab
disp(['Full dynamic simulation computation time:', num2str(timeFullDyn),' s.']);
disp(['Hybrid simulation computation time:', num2str(timeHybrid),' s.']);
```

The complete example can be found in `/example/ex_extended_term_dyn.m`. And the results can also be found in our paper:
* Hybrid QSS and Dynamic Extended-Term Simulation Based on Holomorphic Embedding, arXiv:2104.02877


### Appendix: Variables in settings

#### variable `eventList` 
([back to top](#11-settings-file))
##### Table 1. Definition of `eventList`
Column  | Content
-------| -------------
1 | Event index (can be an integer or a real number)
2 | Event start time
3 | Event end time (no effect for instant event)
4 | Type of event (see [Table 2](#table-2-event-types))
5 | Index of event under its type
6 | Simulation method (default 0.0) (see below [Simulation methods](#simulation-methods))
7 | Timestep (default 0.01)

##### Table 2. Event types
([back to top](#11-settings-file))
Value  | Event type
-------| -------------
0 | Calculate steady-state at start
1 | Add line
2 | Add static load
3 | Add induction motor load
4 | Add synchronous generator
6 | Applying/clearing faults
7 | Cut line
8 | Cut static load
9 | Cut motor load
10| Cut synchronous generator
50| Dynamic process
99| End of simulation

##### Simulation methods
Simulation methods can be specified for each event on the 6th column of `eventList`. It is encoded as a number `x.yz`, where:
* `x` is the method for solving differential equation, where 0 - SAS, 1 - Modified Euler, 2 - R-K 4, 3 - Trapezoidal rule.
* `y` is the method for solving algebraic equation, where 0 - SAS, 1 - Newton-Raphson.
* `z` is whether to use variable time step scheme for numerical integration (`x` is 1, 2 or 3). 0 - Fixed step, 1 - Variable step.

Note that when `x=0`, `y` and `z` are not effective, it automatically uses SAS and variable time steps.

#### variable `bsBus` 
([back to top](#11-settings-file))
Current version only support one black start bus and therefore only the first line will be recognized. Will expand in the future versions.
Column  | Content
-------| -------------
1 | Bus index
2 | Active power of Z component load
3 | Active power of I component load
4 | Active power of P component load
5 | Reactive power of Z component load
6 | Reactive power of I component load
7 | Reactive power of P component load

#### variable `bsSyn`
([back to top](#11-settings-file))
Column  | Content
-------| -------------
1 | Index of synchronous generator
2 | Excitation potential
3 | Active power
4 | Participation factor for power balancing

#### variable `bsInd`
([back to top](#11-settings-file))
Column  | Content
-------| -------------
1 | Index of induction motor
2 | Mechanical load torque

#### variable `Efstd`
([back to top](#11-settings-file))
When there are synchronous generators in the system model, `Efstd` is needed to compute steady state. The `Efstd` is a column vector specifying the excitation potential of every synchronous generator, or it can also be a scalar assigning the excitation potentials of all the generator as the same value.

#### variables `evtLine` and `evtLineSpec`
([back to top](#11-settings-file))
In `eventList`, when the 4th column (event type) equals 1 or 7 (add or cut line, respectively), the index of the line events in `evtLine` corresponds to the 5th column of `eventList`. 
##### variable `evtLine`
([back to top](#11-settings-file))
Column  | Content
-------| -------------
1 | Index of line events (from 5th column of `eventList`)
2 | Start index in `evtLineSpec`
3 | End index in `evtLineSpec`

##### variable `evtLineSpec`
([back to top](#11-settings-file))
Column  | Content
-------| -------------
1 | Index of line 
2 | Add/cut mark, 0 - add line, 1 - cut line
3 | Reserved
4 | Reserved
5 | Reserved

#### variables `evtZip`, `evtZipSpec` and `evtZipSpec2`
([back to top](#11-settings-file))
In `eventList`, when the 4th column (event type) equals 2 or 8 (add/cut static load), the index of the load events in `evtZip` corresponds to the 5th column of `eventList`. 
##### variable `evtLine`
([back to top](#11-settings-file))
Column  | Content
-------| -------------
1 | Index of load events (from 5th column of `eventList`)
2 | Choose `evtZipSpec` (0) or `evtZipSpec2` (1)
3 | Start index in `evtZipSpec` or `evtZipSpec2`
4 | End index in `evtZipSpec` or `evtZipSpec2`

##### variable `evtZipSpec`
([back to top](#11-settings-file))
Column  | Content
-------| -------------
1 | Index of zip loads in system base state 
2 | Add/cut mark, 0 - add load, 1 - cut load

##### variable `evtZipSpec2` (recommended)
([back to top](#11-settings-file))
`evtZipSpec2` has the same format as PSAT ZIP load format, which represents the change of ZIP load. Whether the event is specified as add/cut load does not make difference.

#### variables `evtInd` and `evtIndSpec`
([back to top](#11-settings-file)) 
In `eventList`, when the 4th column (event type) equals 3 or 9 (add or cut induction motors, respectively), the index of the induction motor events in `evtInd` corresponds to the 5th column of `eventList`. 
##### variable `evtInd`
([back to top](#11-settings-file))
Column  | Content
-------| -------------
1 | Index of induction motor events (from 5th column of `eventList`)
2 | Start index in `evtIndSpec`
3 | End index in `evtIndSpec`

##### variable `evtIndSpec`
([back to top](#11-settings-file))
Column  | Content
-------| -------------
1 | Index of induction motor 
2 | Event type, 0 - add motor, 1 - change state, 2 - cut motor
3 | Designated mechanical torque
4 | Designated slip

#### variables `evtSyn` and `evtSynSpec`
([back to top](#11-settings-file))
##### variable `evtSyn`
Column  | Content
-------| -------------
1 | Index of synchronous generator events (from 5th column of `eventList`)
2 | Start index in `evtSynSpec`
3 | End index in `evtSynSpec`

##### variable `evtSynSpec`
([back to top](#11-settings-file))
Column  | Content
-------| -------------
1 | Index of synchronous generator
2 | Event type, 0 - add generator, 1 - cut generator
3 | Designated rotor angle (only effective when adding generator, NaN means the rotor angle is the same with voltage angle).
4 | Designated mechanical power (only effective when adding generator).
5 | Designated excitation potential (only effective when adding generator).

#### variables `evtFault` and `evtFaultSpec`
([back to top](#11-settings-file))
##### variable `evtFault`
Column  | Content
-------| -------------
1 | Index of fault events (from 5th column of `eventList`)
2 | Start index in `evtFaultSpec`
3 | End index in `evtFaultSpec`

##### variable `evtFaultSpec`
([back to top](#11-settings-file))
Column  | Content
-------| -------------
1 | Index of fault line
2 | Position of fault, 0.0 stands for starting terminal and 1.0 stands for ending terminal.
3 | Resistance of fault.
4 | Reactance of fault.
5 | Event type: 0 - add fault; 1 - clear fault.

#### variable `evtDyn`
([back to top](#11-settings-file))
The `evtDyn` variable specifies the indexes of ramping events involving various types of components. 
Column  | Content
-------| -------------
1 | Index of event
2 | Start index in `evtDynPQ`
3 | End index in `evtDynPQ`
4 | Start index in `evtDynPV`
5 | End index in `evtDynPV`
6 | Start index in `evtDynInd`
7 | End index in `evtDynInd`
8 | Start index in `evtDynZip`
9 | End index in `evtDynZip`
10| Start index in `evtDynSh`
11| End index in `evtDynSh`
12| Start index in `evtDynZipRamp`
13| End index in `evtDynZipRamp`
14| Start index in `evtDynTmech`
15| End index in `evtDynTmech`
16| Start index in `evtDynPm`
17| End index in `evtDynPm`
18| Start index in `evtDynEf`
19| End index in `evtDynEf`
20| Start index in `evtDynVref`
21| End index in `evtDynVref`
22| Start index in `evtDynEq1`
23| End index in `evtDynEq1`

#### variable `evtDynPQ` 
([back to top](#11-settings-file))
Column  | Content
-------| -------------
1 | Index of bus
2 | Active power ramping rate
3 | Reactive power ramping rate

#### variable `evtDynPV` 
([back to top](#11-settings-file))
Column  | Content
-------| -------------
1 | Index of bus
2 | Active power ramping rate

#### variable `evtDynInd` 
([back to top](#11-settings-file))
Column  | Content
-------| -------------
1 | Index of induction motor
2 | Mechanical load torque ramping rate

#### variable `evtDynZip` 
([back to top](#11-settings-file))
Column  | Content
-------| -------------
1 | Index of bus
2 | ZIP load ramping rate

#### variable `evtDynSh` 
([back to top](#11-settings-file))
Column  | Content
-------| -------------
1 | Index of bus
2 | Shunt admittance ramping rate

#### variable `evtDynZipRamp` 
([back to top](#11-settings-file))
`evtDynZipRamp` has the same format as PSAT ZIP load format, which represents the ramping direction of ZIP load.

#### variable `evtDynTmech`
([back to top](#11-settings-file)) 
Column  | Content
-------| -------------
1 | Index of synchronous generator
2 | Ramping rate of mechanical power reference value (TG required)

#### variable `evtDynEf` 
([back to top](#11-settings-file))
Column  | Content
-------| -------------
1 | Index of synchronous generator
2 | Ramping rate of excitation potential

#### variable `evtDynVref`
([back to top](#11-settings-file)) 
Column  | Content
-------| -------------
1 | Index of synchronous generator
2 | Ramping rate of exciter reference voltage

#### variable `evtDynEq1` 
([back to top](#11-settings-file))
Column  | Content
-------| -------------
1 | Index of synchronous generator
2 | Ramping rate of transient excitation potential

