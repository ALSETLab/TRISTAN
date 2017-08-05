# TRISTAN: Time seRIes baSed sTability Assessment iNdexes
Methods to assess power system (static and dynamic) stability margins from time-series data

![ScreenShot](https://github.com/ALSETLab/TRISTAN/blob/master/_pics/readmeimage.jpg)

# Introduction and Motivation
The main categories of computations performed for stability analysis purposes in power systems are contingency analysis, and security margin calculations. Contingency analysis methods consist determining system’s response to large disturbance, at a given operating point. Static and dynamic time-domain simulation methods are generally used in contingency analysis to assess the stability of a given power system. Static methods focus on the existence of equilibria and therefore rely on the solution of non-linear algebraic equations that assume equilibrium conditions of the system’s dynamics. Power flow based contingency analysis and continuation power flow (CPF) methods are some examples of static methods. Static methods are usually very efficient but they neither account for post contingency controls that depend on the system’s evolution nor capture more involved instability mechanisms. 

Time-domain methods, on the other hand, may have higher computational demands, but offer higher accuracy and better information (e.g. w.r.t the system’s response to a sequence of events). Such methods are attractive for Dynamic Security Assessment tools, specially those that do not have the facility to exploit the mathematical model's structure and can only use the resulting time-series from time-domain simulations. This is the case of the iPST platform (https://github.com/itesla/ipst) that relies in time-domain simulations for security assessment, and for which the the indexes provided in this repository were first developed.

## Related work
A previoulsy released version of these indexes has been integrated into the iTesla platform, and are available in the open source distribution of the iPST. The available version of these indexes is in the iPST repository: [link](https://github.com/itesla/ipst/tree/e46b47547098915367f4fcfe96301d068b45b2ab/dynamic-indexes).

## Purpose of this repository
Because the indexes are integrated into the platform's workflow, it is difficult to adapt them for other applications and time-series data sources (e.g. simulator outputs different than those in iPST or measurements). Hence, in an effort to provide these basic functionalities for use in other dynamic security assessment platforms or applications, all versions of the indexes developed along with testing data from simulations used in publications is provided in this repository.

## Repository Structure
The repository is organized as follows:
  - ``./_docs/`` includes all the documentation available. The deliverables for the iTesla project for this work, papers written (for which results can be reproduced), presentations, and a synthesis of some of the indexes are provided in sub-folders.
  - ``./_pics/`` provides images with graphical results of the use of the indexes - for illustration purposes.
  - ``./_sour/`` includes all development versions, and final set of indexes. Development versions are in folders market `v0_X...`, and the final version is in a folder named `v1`.
  - ``./_thirdparty`` includes the original code for a [Prony function developed by PNNL](https://github.com/ftuffner/DSIToolbox/blob/master/Ringdown130930/private/prgv2_5.m) which has been extensively modified to fit the purposes of the application. Please refer to this [link](https://github.com/ftuffner/DSIToolbox/) to obtain an updated or newer version of that code.

## Running the different indexes
Each of the indeces are in an independent directories and are implemented in a Matlab function, that in some cases has dependencies to other functions. All functions are tested from individual Matlab scripts that call data and format the required inputs for excecution. All the indexes are located under ``./_sources/``:
  - Static overload index: Run ``static_overload_index_testing.m`` on any of the three data sets ``Over_load-X.mat``.
  - Static overvoltage index: Run ``static_voltage_index_testing.m`` on the data set ``OverUnder_Voltage.mat``.
  - Dynamic Transient index: Run ``dynamic_transcient_index_testing.m`` on the data sets ``transient.mat`` and ``Gentrip_X.mat``.
  - Small signal indexes:
    - Type A: under ``./_source/04_a_ Dynamic Small_Signal_prony_ERA/`` this method allows applying the ERA or Prony method, run ``dynamic_smallsignal_testing.m`` on data set ``SmallSignal.mat``.
    - Type B (Preferred): under ``./_source/04_b_Dynamic_SmallSignal_prony_only/``, a newer version applying only the Prony method that is more stable under unknown inputs, run ``Dynamic_smallsignal_example.m`` on data sets ``7bus_fault.mat`` and, ``KTHNordic32_fault.mat``
  - Voltage stability indes:
    - Type A: under ``./_source/04_b_Dynamic_SmallSignal_prony_only/`` run ``dynamic_voltagestability_testing.m`` on ``LOAD_ALL.mat`` and ``LONNY_ALL.mat``
    - Type B (Preferred): under ``./_source/05_b_Voltage_Stability_Indexes`` run ``Voltage_Stability_index_testing.m`` on data sets ``2bus-system.mat``, ``KTH_Nordic_trip_L32.mat`` and ``KTH_Nordic_trip_L42.mat``.

# Cite our work!
Please cite any of the following three papers depending on the use of our code:

- Static and transient stability indeces:
>> F. R. S. Sevilla and L. Vanfretti, "Static stability indexes for classification of power system time-domain simulations," 2015 IEEE Power & Energy Society Innovative Smart Grid Technologies Conference (ISGT), Washington, DC, 2015, pp. 1-5. [URL:]( http://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=7131846&isnumber=7131775)

- Small signal stability indices:
>> F. R. S. Sevilla and L. Vanfretti, "A small-signal stability index for power system dynamic impact assessment using time-domain simulations," 2014 IEEE PES General Meeting | Conference & Exposition, National Harbor, MD, 2014, pp. 1-5. [URL]( http://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=6938842&isnumber=6938773)

- Voltage Stability Indices:
  - Type A:
>> L. Vanfretti and F. R. S. Sevilla, "A three-layer severity index for power system voltage stability assessment using time-series from dynamic simulations," IEEE PES Innovative Smart Grid Technologies, Europe, Istanbul, 2014, pp. 1-5. [URL]( http://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=7028788&isnumber=7028730).
  - Type B:
>> V.S. Narasimham Arava and L. Vanfretti, "A Method to Estimate Power System Voltage Stability Margins using Time-series from Dynamic Simulations with Sequential Load Perturbations", paper submitted to the IEEE Transactions on Power Systems, 2017. Under Review.

# Developers
Felix Rafael Segundo Sevilla, Venkata Satya Narasimham Arava ([Narasimhamarava](https://github.com/Narasimhamarava)), Luigi Vanfretti ([lvanfretti](https://github.com/lvanfretti))

# License
Thi is free/libre software and the use is completely at your own risk; it can be redistributed and/or modified under the terms of the GNU Public License version 3.

Copyright (C) 2017,  Luigi Vanfretti, Felix Rafael Segundo Sevilla, Venkata Satya Narasimhan Arava.
