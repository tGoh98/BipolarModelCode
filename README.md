# BipolarModelCode
Computational model of Bipolar Cells in Lower Vertebrate Retina. Done for BIOE 401 research (Raphael Lab).

## Files Explained
* ionic_model.mlx
    - Matlab Live Script that computes and plots the voltage and currents.
* calcV.m
    - Computes the voltage membranous voltage.
    - Input:
      - tspan - The time span (s).
      - inj - The amount of current to inject for getI() (pA).
    - Output:
      - t - An array containing the time values for each timestep the ode used.
      - res - A matrix containing the results of the ode for each timestep. Includes intermediate values such as emission  probabilities and Calcium values.
* calcIXY.m
    - Computes the voltage for the respective current.
    - Input:
      - t - Array containing the time values for each timestep the ode used (s).
      - Any intermediate values from calcV() needed to calculate the voltage.
      - V - Array containing the voltage values for each timestep (mV).
    - Output:
      - An array containing the respective current for each timestep (pA).

## Documentation and Comments
In the Matlab console, type ```help [function name]``` to see the documentation for any function.
Constants and parameters are commented with a description and units.

## References
Equations used were based on: 
ADD CITATION HERE
