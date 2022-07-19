# Smagorinsky

The Smagorinsky subgrid-scale (SGS) model was developed by J. Smagorinsky in the
metereological community in 1963[^1]. It is based on the eddy viscosity assumption,
which postulates a linear relationship between the SGS shear stress and the resolved
rate of strain tensor. This model serves as a base for other SGS models.

The subgrid scale stress tensor $\tau_{ij}$ is

$$
\begin{align}
\tau_{ij}&=\overline{u_{i}u_{j}}-\bar{u_{i}}\bar{u_{j}} \\
&=\frac{1}{3}\tau_{kk}\delta_{ij}+\left(\tau_{ij}-\frac{1}{3}\tau_{kk}\delta_{ij}\right)\\
&\approx\frac{1}{3}\tau_{kk}\delta_{ij}-2\nu_{SGS}dev(\overline{D_{ij}})\\
&\approx\frac{2}{3}k_{SGS}\delta_{ij}-2\nu_{SGS}dev(\overline{D_{ij}})
\end{align}
$$

where $\nu_{SGS}$ is the subgrid-scale viscosity and $\overline{D_{ij}}$ is the resolved strain rate tensor 
defined as

$$
\overline{D_{ij}}=\frac{1}{2}\left(\frac{\partial \bar{u_{i}}}{\partial x_{j}}+\frac{\partial \bar{u_{j}}}{\partial x_{i}}\right)
$$

[^1]: Smagorinsky, Joseph. "General circulation experiments with the primitive equations: I. The basic experiment." Monthly weather review 91.3 (1963): 99-164.
