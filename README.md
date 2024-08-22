# CH5140-Process-Modelling-Simulation-and-Analysis
## This project was a part of the course CH5140 Process Modelling Simulation and Analysis

The aim of the project is to model the coupled system of a reactor with a flash downstream considering the effect of transport delay between the two. This problem statement is taken from Pushpavanam et al. (2005) and in this project, the results are reproduced using techniques learnt throughout the course. A detailed presentation and MATLAB codes are uploaded in the repository. The motivation for the work, the methodology used as well the results and conclusions are brielfy presented below.

**Motivation:**

In a typical chemical industry, a reactor is coupled with a seperator downstream to recycle unreacted reactants to the reactor. Hence, these units are coupled to each other by the recycle stream as depicted in Fig 1 which includes a CSTR and a flash unit. We aim to study the stability and dynamics of such a system so as to operate the process at appropriate conditions avoiding unstable states. Further, we consider the effect of transport delay as the stream exiting the reactor does not reach the seperator instantaneously

<div align="center">![Fig 1: The coupled reactor (CSTR) and seperator (flash) system](https://github.com/user-attachments/assets/cdcedf62-fd82-4614-a7b2-32e65503de98)</div>

**Methodology:**
* The governing equations are derived using the species mass balance for each of the unit as well as the entire system and the energy balance for the non isothermal CSTR. To simplify the equations further, we consider the case of a first order irreversible reaction
The further steps are carried out using MATLAB to solve and analyze the system of coupled differential equations with and without delay using the ode23s, dde23 and fsolve functions
* A degree of freedom analysis is performed to identify number of manipulated variables so that their effects on the stability of the system can be analyzed.
* The governing differential equations are non dimensionalized, steady state conditions are calculated and a linear stability analysis is performed for two cases, with and without the transport delay.
* Further bifurcation plots are developed to understand the effect of the manipulated variables to understand their effects on the stability of the steady states.
* The effect of different initial conditions on the dynamics of the system are studied by creating phase plane plots for different values of the non dimensional manipulated variable

**Results:**
* Following is a comparison of the bifurcation plots with Pushpavanam et al. (2005)


* Fig. 4 shows phase plane plots for different values of the non dimensional manipulated variable Da1 depicting its effect on the existence of steady stable states

**Conclusions:**

In case of the system without delay, the instabilities are induced due to the non linearity of the Arrhenius rate law. But, considering the delay we observe new zones of instability solely due to the effect of this lag.
Hence, the delay has an significant effect on the stability and dynamics of thw non isothermal reactor separator system which needs to be considered while designing process plants.


References:

[1] Balasubramanian, P., Pushpavanam, S., Kienle, A., & Balaraman, K. S. (2005). Effect of delay on the stability of a coupled reactor−flash system sustaining an elementary non-isothermal reaction. Industrial & Engineering Chemistry Research, 44(10), 3619-3625. https://doi.org/10.1021/ie040005v 

