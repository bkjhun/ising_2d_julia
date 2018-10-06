# Markov chain Monte Carlo simulation of 2-dimensional ising model
Julia codes for Markov chain Monte Carlo simulation of 2 diemnsional Ising model

## ising_2d_lib.ipynb
Library containing functions used in simulation. Import the library by
```julia
using NBInclude
@nbinclude("ising_2d_lib.ipynb")
```

---------------------------------------

Function that returns 2-dimensional boolean array, random or ordered. Each boolean digit corresponds to up- and down-spin.

**initialize** (Length, initial_state="random")

---------------------------------------

Functions that make spin configuration undergo certain time of stochastic dynamics.

**evolution_rs** (Config, Temperature, Time): 	*Metropolis, random-site flipping*

**evolution_tw** (Config, Temperature, Time):	*Metropolis, typewriter flipping*

Function inputs

>**Config**:	2-dimensional boolean array of spin
>
>**Temperature**:	temperature of the stochastic dynamics
>
>**Time**:	duration of the dynamics. $\text{Time} \times \text{Length}^2$ spin-flipping attempts will be made

---------------------------------------

Functions that return a vector of sampled magnetization and energy data. Each function has distinct Markov chain.

**stream_rs** (Config, Temperature, num_sample, num_sample_burn, period_sample): *Metropolis, random-site flipping*

**stream_tw** (Config, Temperature, num_sample, num_sample_burn, period_sample): *Metropolis, typewriter flipping*

**stream_glauber** (Config, Temperature, num_sample, num_sample_burn, period_sample):   *Glauber, random-site flipping*

**stream_glauber_tw** (Config, Temperature, num_sample, num_sample_burn, period_sample):	*Glauber, typewriter flipping*

Function inputs

>**Config**:	2-dimensional boolean array of spin
>
>**Temperature**:	temperature of the stochastic dynamics
>
>**num_sample**:	number of sample physical properties
>
>**num_sample_burn**:	After this time the stream of sample is taken
>
>**period_sample**:	time between each sample

---------------------------------------

Functions related with autocorrelation.

**autocorrelation** (data, k): 	*k-step autocorrelation of data*

**autocorrelation_function** (data, kmax):	*Returns 1- to kmax-step autocorrelation of data*

**autocorrelation_time** (data, kmax):	*Returns autocorrelation time, which is the first time that autocorrelation is less than exp(-1)*

## data_analysis.ipynb
This notebook uses previously computed data of correlation time, average and variance of magnetization, average and variance of energy.

## ising_2d_timing.ipynb
This notebook measures time required to perform MCMC.

# Results
Random-site flipping dynamics has around x4 correlation time compared to typewriter flipping dynamics.

## Critical exponents

* &nu; = 1.0
* &beta; = 0.125
* &gamma; = 1.75
