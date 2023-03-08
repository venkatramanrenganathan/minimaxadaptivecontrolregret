# Regret Analysis for Minimax Adaptive Controller
This repository contains the MATLAB code for studying the regret of a minimax adaptive control algorithm. 

**Associated Paper:** Venkatraman Renganathan, Andrea Ianelli, and Anders Rantzer, `Online Learning and Regret Analysis for Minimax Adaptive Control Algorithm`, Submitted to IEEE Conference on Decision and Control, Singapore, 2023.

# Dependencies
- Matlab
- Yalmip Toolbox
- MOSEK solver

# Procedure to run the code
1. Run the matlab code `RegretDoubleIntegrator.m` which will load the required system and controller data and generate the desired plots.

## Variations while running `RegretDoubleIntegrator.m` file
    * Set `disturbanceSelect = 1` in line `137` of `RegretDoubleIntegrator.m` for simulations with Hinfinity adversarial disturbance
    * Set `disturbanceSelect = 2` in line `137` of `RegretDoubleIntegrator.m` for simulations with confusing adversarial disturbance
    * Set `disturbanceSelect = 3` in line `137` of `RegretDoubleIntegrator.m` for simulations with sinusoidal adversarial disturbance
    * Set `addrandomFlag = 1` in line `139` of `RegretDoubleIntegrator.m` for simulations with disturbance with randomness too.

# Contributing Authors
1. [Venkatraman Renganathan - Lund University](https://github.com/venkatramanrenganathan)
2. [Andrea Ianelli - University of Stuttgart](https://andreaian.github.io)
3. [Anders Rantzer - Lund University](https://control.lth.se/personnel/personnel/anders-rantzer/)

# Funding Acknowledgement
V. Renganathan and A. Rantzer have been supported by the *European Research Council (ERC)* under the European Union’s Horizon 2020 research and innovation program under grant agreement No `834142 - Scalable Control`.

# Affiliation
1. [Department of Automatic Control, Lund University, Sweden](https://control.lth.se)
2. [Institute for Systems Theory and Automatic Control (IST) at University of Stuttgart, Germany](https://www.ist.uni-stuttgart.de)
