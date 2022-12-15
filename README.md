# Shallow water equations with Multigrid method

This model solving a 2D shallow water equations refering to Spitaleri and Corinaldesi(1997). Semi-implicit finite difference method is used for the time discretizing, and Semi-Lagrange method is used for the convection term discretization. When solving a Poisson-like equations, so-called "Multigrid method" is used. Test problem is about the behaviour of basin water under time variance tide.

## Equations
$$
\frac{\partial u}{\partial t} + u \frac{\partial u}{\partial x} + v \frac{\partial u}{\partial y} + g \frac{\partial z}{\partial x} = - \gamma u + f v
$$

$$
\frac{\partial v}{\partial t} + u \frac{\partial v}{\partial x} + v \frac{\partial v}{\partial y} + g \frac{\partial z}{\partial y} = - \gamma v - f u
$$

$$
\frac{\partial z}{\partial t} + \frac{\partial [(h+z)u]}{\partial x} + \frac{\partial [(h+z)v]}{\partial y}  = 0
$$

## Some notations for Multigrid methods
(a) coarser-grid construction: standard coarsening

(b) type of grids: Arakawa C-grid

(c) relaxation for error smoothing: Gauss-Seidel Lexicographic order

(d) restriction: local averaging

(e) interpolation: weighted interpolation

(f) multigrid cycling: V-cycle