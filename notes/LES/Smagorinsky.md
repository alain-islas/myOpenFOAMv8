# Smagorinsky

The Smagorinsky subgrid-scale (SGS) model was developed by J. Smagorinsky in the
metereological community in 1963. It is based on the eddy viscosity assumption,
which postulates a linear relationship between the SGS shear stress and the resolved
rate of strain tensor. This model serves as a base for other SGS models.

\usepackage{amsmath}

The subgrid scale stress tensor $\tau_{ij}$ is

$$
\begin{align}
\tau_{ij}&=\overline{u_{i}u_{j}}-\bar{u_{i}}\bar{u_{j}} \\
&=\frac{1}{3}\tau_{kk}\delta_{ij}+\left(\tau_{ij}-\frac{1}{3}\tau_{kk}\delta_{ij}\right)\\
&\approx\frac{1}{3}\tau_{kk}\delta_{ij}-2\nu_{SGS}dev(\overline{D_{ij}})
\end{align}
$$
