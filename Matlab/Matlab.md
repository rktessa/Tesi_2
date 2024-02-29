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
- sirs_0 working and simulate a sisrs. With the coefficients $\beta = \gamma * 2 * N$ and $\gamma = 0.2679$ and $\delta = 1/205$ the peaks of infected are temporaly similar to the covid in italy durign 2020. 
- sirs1 try to correct the ode not working simulaiton, but firt implementing in runge kutta a variable beta.  