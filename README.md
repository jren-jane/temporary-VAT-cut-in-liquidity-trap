# temporary-VAT-cut-in-liquidity-trap
This repository performs a perfect foresight simulation of the temporary VAT cut in the UK in 2009.

## Introduction to the model
The tax cut is meaningful to study due to its temporary nature and the liquidity trap. As central banks are running out of monetary measures, the expectations for the fall and rise in inflation could possibly lower the real interest rate through the Fisher equation, and stimulate the economy.

The code, written in dynare, is based on a large-scale dynamic stochastic general equilibirum model for which the parameters are obtained from Bayesian estimation. It generates impulse responses for critical macroeconomic variables such as output, inflation, unemployment rate, consumption, investment, etc.  

## File description
### Estimation
- Posterior_Generation.mod is for estimating the model parameters by Bayesian methods.
- data.mat is processed data input of estimation

### Deterministic simulation
All the other files are for deterministic simulation. 

## Implementation
You need to install dynare to run the code.

## Results
![img](https://github.com/jren-jane/temporary-VAT-cut-in-liquidity-trap/blob/master/output.png)
