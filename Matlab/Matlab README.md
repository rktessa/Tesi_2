# Notes on Matlab codes

## List of files
### Update at 14.11.23
- sir_1.m contain simulation  of a simple SIR model while varying the beta and gamma coefficients
- opinion_0.m it implements a Careless, Compliant, Anti model for the behaviour of agent after acquiring information about the disease

## Update at 21.11.23
 Creation of a class of functions "functionsContainer.m" to call the created function from all files. It contains the function created in sir_1.m and opinion_0.m, but can be used with N = 1, instead of N = 10000. 

## Update 29.02.24
- elragiandonofrio dimple function erlagian to simulate queueing ?
- sis_0 implented to try simulate a sis
- sirs_0 working and simulate a sisrs. With the coefficients $\beta = \gamma * 2 * N$ and $\gamma = 0.2679$ and $\delta = 1/205$ 
 the peaks of infected are temporaly similar to the covid in italy durign 2020. 
- sirs1 try to correct the ode not working simulation, but first implementing in runge kutta a variable beta.  
- sirs_multiple_sensitivity here, multiple simualations are performed to see how vary the situation at coefficient transformation. 

## Update 24.03
- simulated multiple behavioural models, with several simulations done varying the coefficients
- analysis of the result and work on the dedicated chapter
- epi_behavioural_model_00 file contain the multi-layer network

## Update 29.04
- There are several new files, in particular the analysis of behavioral and sirs model, that are done with live script for improve the comprehension. 
- The folder is been completely re-organised, all the files are now divided in several subfolders depending on their topic.
- With this modification the files that use the "functionsContainer" must have a new line of code for identify the position of the file.
- the command: 
    addpath('../Matlab')
    solve the problem

## Update 12.12
- In the Multilayer system folder there are all the code to reproduce the calculations and pictures realized in the thesis.