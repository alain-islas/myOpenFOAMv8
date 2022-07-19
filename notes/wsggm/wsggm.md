# Weighted-sum of gray gas model (WSGGM)

The weighted sum of gray gas model (WSGGM) was first developed by
Hottel and Sarofim[^1]. It replaces the spectrum with few gray gases
and a transparent window according to

$$
\varepsilon = \sum_{i=0}^{N_{g}}a_{\varepsilon,i}\left(T\right)\left[1-\exp{(-\kappa_{i}p_{a}L)}\right]
$$

where $N_{g}$ is the number of gray gases and $a_{i}$ are the emissivity weighting factors.
The bracketed quantity in Eq. () is the $i$-th gray gas emissivity with
banded absorption coefficient $\kappa_{i}$ and pressure-path length $p_{a}L$. The pressure $p_{a}$
is expressed by summing the partial pressures of the participating gases, namely $H_{2}O$ and $CO_{2}$

$$
p_{a}=(X_{CO_{2}}+X_{H_{2}O})p
$$

where $X_{i}$ denotes the molar fraction of each species and $p$ is the total pressure in atm.
To represent the transparent parts of the spectrum, the banded absorption coefficient $\kappa_{i=0}=0$.
Since total emissivity approaches unity in the limit of the pressure-path length, the emissivity
weighting factors must sum unity, and all adopt positive values. This implies that 
$a_{\varepsilon,0}=1-\sum a_{\varepsilon,i}$ such that only $N_{g}$ weighting
factors need to be determined.

Commonly, the emissivity weighting factors are assumed to be a temperature dependent polynomial
function of order $(N_{g}-1)$, i.e.

$$
a_{\varepsilon,i}\left(T\right)=\sum_{j=1}^{N_{g}}b_{\varepsilon,i,j}T^{j-1}
$$

where $b_{\varepsilon,i,j}$ are the polynomial coefficients. In the standard form of (Smith et al., 1982)[^2]
these coefficients are valid for specific molar ratios, $MR = X_{H_{2}O}/X_{CO_{2}}=1.0\,2.0$ and preferably should
be used only if the composition of the mixture is known.

[^1]: Hottel, Hoyt C. "Radiant heat transmission." WH McAdams. Heat Transmission (1954).
[^2]: Smith, T. F., Z. F. Shen, and J. N. Friedman. "Evaluation of coefficients for the weighted sum of gray gases model." (1982): 602-608.
