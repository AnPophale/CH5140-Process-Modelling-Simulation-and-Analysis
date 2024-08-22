# CH5140-Process-Modelling-Simulation-and-Analysis
### Modelling of a coupled reactor-separator system with transport delay

The aim of the project is to model the coupled system of a reactor (CSTR) with a separator (Flash) downstream considering the effect of transport delay between the two. This problem statement is taken from Pushpavanam et al. (2005), and in this project the results are reproduced using techniques learnt throughout the course. A [detailed presentation](https://github.com/AnPophale/CH5140-Process-Modelling-Simulation-and-Analysis/blob/main/PMSA%20Term%20Paper%20Presentation.pdf) and [MATLAB codes](https://github.com/AnPophale/CH5140-Process-Modelling-Simulation-and-Analysis/tree/8ffd7eba8713dd2d7d0d88d01121c08533d17a86/MATLAB%20Code) are uploaded to this repository. The motivation for the work, the methodology used as well as the results and conclusions are briefly presented below.

**Motivation:**  
In a typical chemical industry, a reactor is coupled with a separator downstream to recycle unreacted reactants to the reactor. Hence, these units are coupled to each other by the recycle stream, as depicted in Fig. 1 which includes a CSTR and a flash unit. We aim to study the stability and dynamics of such a system so as to operate the process at appropriate conditions avoiding unstable states. Further, we consider the effect of transport delay as the stream exiting the reactor does not reach the separator instantaneously.  
<p align="center">
  <img src="https://github.com/user-attachments/assets/cdcedf62-fd82-4614-a7b2-32e65503de98" alt="Fig 1: The coupled reactor (CSTR) and separator (Flash) system" style="width: 50%;">
</p>
<p align="center">
  <em>Figure 1: The coupled reactor (CSTR) and separator (Flash) system</em>
</p>


**Methodology:**
* The governing equations are derived using the species mass balance for each of the unit as well as the entire system and the energy balance for the non-isothermal CSTR. To simplify the equations further, we consider the case of a first order irreversible reaction.
  
The further steps are carried out using MATLAB to solve and analyze the system of coupled differential equations with and without delay using the ode23s, dde23 and fsolve functions.
* A degree of freedom analysis is performed to identify number of manipulated variables so that their effects on the stability of the system can be analyzed.
* The governing differential equations are non dimensionalized, steady state conditions are calculated, and a linear stability analysis is performed for two cases, with and without the transport delay.
* Further bifurcation plots are developed to understand the effect of the manipulated variables to understand their effects on the stability of the steady states.
* The effect of different initial conditions on the dynamics of the system are studied by creating phase plane plots for different values of the non dimensional manipulated variable.

**Results:**
* Following is a comparison of the bifurcation plots with Pushpavanam et al. (2005)  
<p align="center">
  <img src="https://github.com/user-attachments/assets/24e76a52-e786-420b-895e-af25e9506e35" alt="Fig 2: Bifurcation plots for the system without delay" style="width: 75%;">
</p>
<p align="center">
  <em>Figure 2: Bifurcation plots for the system without delay</em>
</p>

<p align="center">
  <img src="https://github.com/user-attachments/assets/202c2f49-45f9-497d-aaa4-494607e58fe4" alt="Fig 3: Bifurcation plots for the system with delay" style="width: 75%;">
</p>
<p align="center">
  <em>Figure 3: Bifurcation plots for the system with delay</em>
</p>

* Fig. 4 shows phase plane plots for different values of the non dimensional manipulated variable Da1 depicting its effect on the existence of steady stable states  
<p align="center">
    <img src="https://github.com/user-attachments/assets/8dfa51e1-3c09-4e9c-8055-2e3773345277" style="width: 75%;">
</p>
<p align="center">
  <em>Figure 4: Phase plane plots without delay for different values of the non-dimensional manipulated variable Da1</em>
</p>


**Conclusions:**  
In the case of the system without delay, the instabilities are induced due to the non linearity of the Arrhenius rate law. But, considering the delay we observe new zones of instability solely due to the effect of this lag.
Hence, the delay has a significant effect on the stability and dynamics of the non isothermal reactor separator system which needs to be considered while designing process plants.


References:  
[1] Balasubramanian, P., Pushpavanam, S., Kienle, A., & Balaraman, K. S. (2005). Effect of delay on the stability of a coupled reactor−flash system sustaining an elementary non-isothermal reaction. Industrial & Engineering Chemistry Research, 44(10), 3619-3625. https://doi.org/10.1021/ie040005v 

