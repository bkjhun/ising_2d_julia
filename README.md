# ising_2d_julia
Julia codes for Markov chain Monte Carlo simulation of 2 diemnsional Ising model

* **data_analysis.ipynb** uses previously computed data of correlation time, average and variance of magnetization, average and variance of energy.
* **ising_2d_timing.ipynb** measures time required to perform MCMC.

# Theory

\begin{equation}
Z = \sum e^{-\beta E[\{\sigma_i\}]} = \sum e^{-\beta (J \sum_{<i,j>} \sigma_i \sigma_j - h \sum_{i} \sigma_i)}
\end{equation}

\begin{equation}
-\frac{\partial}{\partial\beta} \ln Z = \sum E[\{\sigma_i\}] \frac{1}{Z} e^{-\beta E[\{\sigma_i\}]} = \left<E \right>
\end{equation}

\begin{align}
c_V &= \frac{1}{N} \frac{\partial}{\partial T} \left<E \right> = \frac{\partial \beta}{\partial T} \left( - \frac{\partial^2}{\partial\beta^2} \ln Z \right) \\
    &= \frac{1}{N T^2} \sum E[\{\sigma_i\}]^2 \frac{1}{Z} e^{-\beta E[\{\sigma_i\}]} - \frac{1}{N T^2} \sum E[\{\sigma_i\}] \frac{1}{Z} e^{-\beta E[\{\sigma_i\}]} \sum E[\{\sigma_i\}] \frac{1}{Z} e^{-\beta E[\{\sigma_i\}]} \\
    &= \frac{1}{N T^2} \left( \left<E^2 \right> - \left<E \right>^2 \right) = \frac{N}{T^2} \left( \left<\left(E/N \right)^2 \right> - \left<E/N \right>^2 \right)
\end{align}

\begin{equation}
\frac{\partial}{\partial h} \ln Z = \sum \beta \left( \sum_{i} \sigma_i \right) \frac{1}{Z} e^{-\beta (J \sum_{<i,j>} \sigma_i \sigma_j -h  \sum_{i} \sigma_i)} = \frac{1}{T} \left<M \right>
\end{equation}

\begin{align}
\chi &= \frac{1}{N} \frac{\partial}{\partial h} \left<M \right> = \frac{\partial}{\partial h} \left( T \frac{\partial}{\partial h} \ln Z \right) \\
    &= \frac{1}{NT} \sum \beta^2 \left( \sum_{i} \sigma_i \right)^2 \frac{1}{Z} e^{-\beta (J \sum_{<i,j>} \sigma_i \sigma_j -h \sum_{i} \sigma_i)} - \frac{1}{NT} \sum \beta \left( \sum_{i} \sigma_i \right) \frac{1}{Z} e^{-\beta (J \sum_{<i,j>} \sigma_i \sigma_j -h \sum_{i} \sigma_i)}  \sum \beta \left( \sum_{i} \sigma_i \right) \frac{1}{Z} e^{-\beta (J \sum_{<i,j>} \sigma_i \sigma_j -h \sum_{i} \sigma_i)} \\
    &= \frac{1}{NT} \left( \left< M^2 \right> - \left< M \right>^2 \right) = \frac{N}{T} \left( \left< m^2 \right> - \left< m \right>^2 \right)
\end{align}s

# Results
Random-site flipping dynamics has around x4 correlation time compared to typewriter flipping dynamics.

## Critical exponents

* &nu; = 1.0
* &beta; = 0.125
* &gamma; = 1.75
