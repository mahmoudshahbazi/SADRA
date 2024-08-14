# SADRA
## A Universal AC/DC Branch Model for Optimal Power Flow Studies

SADRA is a new efficient universal AC/DC branch model for modelling hybrid AC/DC systems for optimal power flow studies with provisions of voltage and power controls. It  provides a framework for modelling wide variety of AC, DC and AC/DC elements including VSCs and VSC-interfaced elements (including point-to-point and multi terminal HVDC), phase-shifter and tap-changing transformers, and in general, hybrid AC/DC systems, in one compact model. SADRA provides a direct link between AC and DC grids, and therefore it is able to use conventional AC equations for modelling. Moreover, it is capable of implementing VSC control actions as well. Due to its compact and simple structure, SADRA is fast and robust.


SADRA is implemented in AIMMS, which is is an optimisation software that supports a wide range of mathematical optimization problems and provides access to multiple solvers such as IPOPT and CONOPT. AIMMS is a powerful tool for optimisation and its user friendly graphical interface makes it easy to use for implementation of optimisation problems, with minimum coding required. 
 
 
## How to Start
 
A fully functional model of AC OPF and its explanation is previosuly published in AIMMS Academy website and is available for the interested readers in [this link](https://how-to.aimms.com/Articles/510/opf.html).

## The provided model
To prove SADRA's performance, versatility and speed, a large AC/DC systems with 3120 AC and 5 DC buses is implemented here. 
