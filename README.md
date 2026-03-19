## Building & Running

`make` builds the code.

`./build/exe` runs the simulation, dumping all output images into the current directory.

## Equations of motion

From [A literate lattice Boltzmann code](http://literatelb.org/#org0ae0da8):

We evolve a distribution function $f(x, v, t)$ where $v = dx/dt$ with governing eqn where the velocity is discretised into (in 2D) 9 components, $\xi_i$. 

$$
(\partial_t + v \cdot \nabla_x)f = \Omega(f),
$$
where $\Omega(f)$ is a collision operator, usually taken to be the BGK operator,
$$
\Omega(f) = - \frac{f - f^{eq}}{\tau},
$$
where $tau$ is some relaxation parameters and
$$
f_i^{eq}(x, t) = \omega_i \rho (1 + u \cdot \xi_i / c_s^2 + (u \cdot \xi_i)^2 / (2 c_s^4) - |u|^2 / (2 c_s^2))
$$

The macroscopic variables are given by moments over $\xi_i$:

$$
\rho = \sum_{k=0}^{q-1}f_k;\n
u = (\sum_{k=0}^{q-1} \xi_k f_k)/\rho
$$
