# Regret Analysis for Minimax Adaptive Controller
This repository contains the MATLAB code for studying the regret of a minimax adaptive control algorithm. 

**Associated Paper:** Venkatraman Renganathan, Andrea Iannelli, and Anders Rantzer, `Online Learning and Regret Analysis for Minimax Adaptive Control Algorithm`, Submitted to IEEE Conference on Decision and Control, Singapore, 2023.

# Dependencies
- Matlab
- Yalmip Toolbox
- MOSEK solver

# Procedure to run the code
1. Run the matlab code `RegretRandomSystemDemo.m` which will load the required system and controller data and generate the desired plots.

## Variations while running `RegretDoubleIntegrator.m` file
    * Set `disturbanceSelect = 1` in line `175` for simulations with Hinfinity adversarial disturbance
    * Set `disturbanceSelect = 2` in line `175` for simulations with confusing adversarial disturbance
    * Set `disturbanceSelect = 3` in line `175` for simulations with sinusoidal adversarial disturbance

## Compute or load the pre-computed gains of minimax adaptive control and Hinfinity control
   * Set `computeFlag = 1` in line `54` to compute minimax adaptive control and Hinfinity control gains
   * Set `computeFlag = 0` in line `54` to load the pre-computed minimax adaptive control and Hinfinity control gains
   
## Set the index for true system model which will generate the data
   * Set `modelNum` to a value in the range `1`-`4` in line `178` 
   
## Set the index for system model which confusing adversarial disturbance should trick the policy to choose
   * Set `cheatModelNum` to a value in the range `1`-`4` in line `181` 
   
## Intricacies While Simulating Regret
   * When `disturbanceSelect = 1` or `disturbanceSelect = 3` in line `175`, simulate Hinfinity system first and then use the disturbance generated for the minimax adaptive control system
   * When `disturbanceSelect = 2` in line `175`, simulate the minimax adaptive control system first and then use the disturbance generated for the Hinfinity system


# Contributing Authors
1. [Venkatraman Renganathan - Lund University](https://github.com/venkatramanrenganathan)
2. [Andrea Iannelli - University of Stuttgart](https://andreaian.github.io)
3. [Anders Rantzer - Lund University](https://control.lth.se/personnel/personnel/anders-rantzer/)

# Funding Acknowledgement
V. Renganathan and A. Rantzer have been supported by the *European Research Council (ERC)* under the European Union’s Horizon 2020 research and innovation program under grant agreement No `834142 - Scalable Control`.

# Affiliation
1. [Department of Automatic Control, Lund University, Sweden](https://control.lth.se)
2. [Institute for Systems Theory and Automatic Control (IST) at University of Stuttgart, Germany](https://www.ist.uni-stuttgart.de)
