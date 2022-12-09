# Shallow water equations with Multigrid method

This model solving a 2D shallow water equations refering to Spitaleri and Corinaldesi(1997). Semi-implicit finite difference method is used for the time discretizing, and Semi-Lagrange method is used for the convection term discretization. When solving a Poisson-like equations, so-called "Multigrid method" is used. Test problem is about the behaviour of basin water under time variance tide.

# Equations
$$
\frac{\partial u}{\partial t} + u
$$
