# PowerSAS.m

**PowerSAS.m** is a robust, efficient and scalable power grid analysis framework based on semi-analytical solutions (SAS) technology. The **PowerSAS.m** is the version for MATLAB/Octave users. It currently provides the following functionalities (more coming soon!):

* **Steady-state analysis**, including power flow (PF), continuation power flow (CPF), contingency analysis.
* **Dynamic security analysis**, including voltage stability analysis, transient stability analysis, and flexible user-defined simulation.
* **Hybrid extended-term simulation** provides adaptive QSS-dynamic hybrid simulation in extended term with high accuracy and efficiency.

### Key features
* **High numerical robustness.** Backed by the SAS approach, the PowerSAS tool provides much better convergence than the tools using traditional Newton-type algebraic equation solvers when solving algebraic equations (AE)/ordinary differential equations (ODE)/differential-algebraic equations(DAE).  
* **Enhanced computational performance.** Due to the analytical nature, PowerSAS provides model-adaptive high-accuracy approximation, which brings significantly extended effective range and much larger steps for steady-state/dynamic analysis. PowerSAS has been used to solve large-scale system cases with 200,000+ buses.
* **Customizable and extensible.** PowerSAS supports flexible customization of grid analysis scenarios, including complex event sequences in extended simulation term. 

### Citing

If you are using **PowerSAS.m** in research work to be published, please include explicit citation of our work in your publication. Please place the following entries in your bibliography:

[1]. J. Liu, R. Yao, F. Qiu, Y. Liu and K. Sun, "PowerSAS.m – An Open-Source Power System Simulation Toolbox Based on Semi-Analytical Solution Technologies," in IEEE Open Access Journal of Power and Energy, In Press.

The corresponding BiBTex citations are given below:

    @ARTICLE{powersas,
             author={Liu, Jianzhe and Yao, Rui and Qiu, Feng and Liu, Yang and Sun, Kai},
             journal={IEEE Open Access Journal of Power and Energy}, 
             title={PowerSAS.m – An Open-Source Power System Simulation Toolbox Based on Semi-Analytical Solution Technologies}, 
             year={2023},
             doi={10.1109/OAJPE.2023.3245040}
            }
