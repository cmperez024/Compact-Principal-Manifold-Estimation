# Fixed-design PME and extrinsic local regression

This experiment compares both estimators in Lin et al. (2017): the extrinsic
kernel estimator in Equations (3)-(4) and the extrinsic local-linear estimator
in Equations (5)-(7). Both calculate an ambient Euclidean regression estimate
and then project it onto the PME manifold.

The fixed-design data-generating model is

\[
t_i \stackrel{\mathrm{once}}{\sim} \operatorname{Unif}(-1,1),\qquad
u_i = \frac{1+\sin(\pi t_i)}{2}\pmod 1,
\]

\[
\mu_i=F(u_i),\qquad
X_i=\mu_i+\varepsilon_i,\qquad
\varepsilon_i\overset{\mathrm{iid}}{\sim}N_2(0,0.08^2I_2),
\]

where

\[
F(u)=\{1+0.30\cos(10\pi u)\}
      \begin{pmatrix}\cos(2\pi u)\\ \sin(2\pi u)\end{pmatrix}.
\]

After generating the design, all analysis conditions on the realized values
of `t_i`. There is no stochastic latent phase and hence no measurement-error
model for the predictor or projection index. The nonmonotone deterministic
map from `t_i` to `u_i` also means the observed predictor is not itself the
flower's projection index.

Five-fold validation jointly selects the PME smoothing parameter and regression
bandwidth separately for the kernel and local-linear estimators. Each validation
response remains raw and is not PME-denoised. PME fits and the complete final
RDS are cached so plotting changes do not require refitting.

`plot-selected-pme.R` recreates the figure from the saved RDS. Its second and
third panels show the complete PME manifolds selected by the kernel and
local-linear downstream criteria, rather than their regression predictions.

## Result

Both estimators selected `lambda = 1e-8` and `h = 0.01` from grids
`1e-14,...,1e-4` and `0.004,...,0.38`. All 55 fold-specific PME fits
converged. Mean five-fold out-of-fold RMSE against the untouched noisy
responses was `0.12420` for the extrinsic kernel estimator and `0.12015` for
the extrinsic local-linear estimator. The corresponding simulation-oracle
RMSE against the conditional mean was `0.05508` and `0.04736`, respectively.

The raw observed responses are not manifold-valued because their error is
ambient. Lin et al.'s estimators are applied only after PME projects the
training responses onto the estimated manifold. This is therefore the
advisor's proposed PME-preprocessing framework, rather than a direct
simulation under Lin et al.'s assumption that the observed responses already
lie on a known manifold.
