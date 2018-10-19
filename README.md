Markov chain Monte Carlo simulation and finite-size scaling of 2 diemnsional Ising model in Julia language.

## ising_2d.ipynb
* visualizes of spin configuration of temperature corresponding to ferromagnetic, paramagnetic, and critical.
* plots autocorrelation function of four MCMC methods: (random-site flipping, typewriter flipping) &bigotimes; (Metropolis dynamics, Glauber dynamics)
* calculates correlation time of systems of various size and temperature, and saves the data
* calculates mean and variance of energy and magnetization, and saves the data

### Computation time
Tested on:  
Intel(R) Core(TM) i7-3770 CPU @ 3.40GHz 3.40GHz

* Visualization: **0.3** s
* Autocorrelation function plot: **9** s
* Autocorrelation function plot of various methods: **259** s = 4 min 19 s
* Correlation time calculation (Typewriter): **4579** s = 1 h 16 min
* Correlation time calculation (Random-site): **22341** s = 6 h 12 min
* Physical quantity calculation: **8210** s = 2 h 17 min

Correlation time of random-site flipping dynamics is not required unless you want to investigate correlation times of various MCMC dynamics.

## data_analysis.ipynb
plots graphs using data previously computed from **ising_2d.ipynb**
* correlation time
* magnetization, magnetic susceptibility, energy, specific heat of systems with various size and temperature
* finite-size scaling of magnetization and magnetic susceptibility

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
>**num_sample_burn**:	After this time, samples are taken. Must be taken long enough so that the system is thermalized after this time.
>
>**period_sample**:	time between each sample

---------------------------------------

Functions related with autocorrelation.

**autocorrelation** (data, k): 	*k-step autocorrelation of data*

**autocorrelation_function** (data, kmax):	*Returns 1- to kmax-step autocorrelation of data*

**autocorrelation_time** (data, kmax):	*Returns autocorrelation time, which is the first time that autocorrelation is less than exp(-1)*

## ising_2d_timing.ipynb
This notebook measures time required to perform MCMC.

# Results
Random-site flipping dynamics has around 420% correlation time compared to typewriter flipping dynamics.

## Critical exponents

* &nu; = 1.0
* &beta; = 0.125
* &gamma; = 1.75
