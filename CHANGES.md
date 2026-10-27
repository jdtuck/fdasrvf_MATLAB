# 3.7.X
* `DP1` alignment (`DynamicProgrammingQ`, also used for curves) resamples q1 and q2 three times finer for its cost (`scl = 3`, previously 1, settable as an optional fifth argument) and interpolates q2 linearly instead of rounding to the nearest point; warps are less jagged and fit better (on a known warp: within 0.26 instead of 0.63 grid cells, alignment error 0.85x instead of 1.35x the true warp's), at about 2.8x the run time
* less noisy warping maps and phase distances: `gam_to_psi` (and so `gam_to_v`) rescales `psi` to unit norm so it stays on the Hilbert sphere after smoothing/clipping; `inv_exp_map`, `elastic_distance` and `elastic_distance_curve` compute the angle with `atan2` instead of `acos`, which was only accurate to ~1e-8 near zero. `gam_to_psi`, `gam_to_v` and `gam_to_h` now default to `smooth=false`: the fixed `SmoothingParam` smoothed more or less depending on the grid size, flattened DP warps (understating phase distances), did not invert through `psi_to_gam`/`v_to_gam`/`h_to_gam`, and needed the Curve Fitting Toolbox. This changes `fdajpca`/`fdahpca`, which use the default; `gam_to_h(smooth=false)` now clamps `gamdot` at `eps` before the `log`, as the smoothed path already did
* `f_to_srvf` with `smooth=false` now differentiates the interpolating cubic spline (the `spl` argument is ignored), and `srvf_to_f` integrates the spline of `q|q|` instead of using `cumtrapz`; the round trip `f -> q -> f` is now O(h^4) accurate (previously O(h^2) and poor for oscillatory functions). `smooth` now defaults to `false` package-wide (the interpolating spline); the alignment classes (`fdawarp`, `pcaTB`, `bootTB`, ...) honor `option.spl` (default `true`) and use the smoothing spline only when `option.spl = false`. `smooth=true` remains lossy by design
* fix multinomial logistic class ordering: training and prediction now both use the standard softmax (`elastic_mlpcr_regression`, `elastic_mlogistic`), and remove the undefined `phi` call
* fix `elastic_mlpcr_regression` crashes (undefined `N`/`N1`, overwritten `alpha`/`b`, nonexistent `LL`/`y_labels` properties, undefined `y_labels`/`option.y` in accuracy)
* `align_fPCA` applies the final warping before centering so `fn`/`qn` match the returned `gam`
* clip inner products to [-1,1] before every `acos`
* clamp warping derivatives at zero before taking `sqrt`
* guard zero-velocity samples in `findBasisNormal` and `Basis_Normal_A`; `findBasisNormal`/`projectTangent` now support curves in R^n and the normal basis is exactly orthogonal to q under `InnerProd_Q`
* fix `Gram_Schmidt`, `Project_Tangent` and `Parallel_Transport_C` (previously errored)
* guard zero distances in `SqrtMedian` and `time_warping_median`
* `Find_Best_Rotation` chooses the reflection correction from `det(U*V')`
* `predict(newdata)` in `elastic_pcr_regression`, `elastic_lpcr_regression` and `elastic_mlpcr_regression` now projects new data with the fPCA class's `project` method (previous inline projection errored); scores use matrix products instead of broadcasting row by column
* fix `project` in `fdajpca`, `fdavpca`, `fdahpca` (`optimum_reparam` argument order, `gradietn` typo, alignment now uses the training `lambda`/`method`)
* `fdahpca` coefficients have one row per function (were sized by the number of time points)
* `fdajpca.calc_fpca` honors the `id` argument instead of always using the midpoint
* fix `fdahpns`: `calc_fpns` passed `stds` instead of the PNS fit to `fastPNSe2s`, and `project` had the `optimum_reparam` argument order wrong, ignored the training `lambda`/`method`/`penalty`, and used the wrong PNS map (now `fastPNSs2e`)
* `elastic_lpcr_regression` constructor was misnamed; fix undefined `N1`, `LL`, `y_labels` in the PCR regression classes; pass alignment options to `time_warping` as name-value arguments
* add `RegressionPredictTest`
* `time_warping` and `time_warping_median`: `fn`/`qn` are now the original data warped by the returned `gam` (previously they kept the preceding iteration's alignment while `gam` came from the final matching step); `mqn` is the mean of `qn` for `time_warping`
* `time_warping_median` honors `smooth` (the SRVFs were computed from the unsmoothed data)
* `fdahpca.calc_fpca` with `log_der = true` returned all-NaN `gam_pca` (a row vector was passed to `h_to_gam`); the direction buffer is now sized by `stds` and `no`
* add `TimeWarpingConsistencyTest`
* `optimum_reparam` defaults to the `roughness` penalty for every method, including `DP1`, so `lambda > 0` changes the default alignment again (the `none` default for `DP1` from 3.7.0 made `time_warping(lambda)` ignore `lambda` once the MEX files were rebuilt). For `DP`/`DP1` this is the integral of (1-sqrt(gamdot))^2, the same default as R (`roughness`) and Python (`l2psi`)
* `time_warping`, `time_warping_median`, `multiple_align_functions` and `ppd` take `option.penalty` (default `roughness`) and store it in the new `penalty` property; `pairwise_align` takes a `penalty` argument
* `DynamicProgrammingQ` (`DP1`) is back to the original grid (`scl = 1`); 3.7.0 accidentally reverted it to the 5x upsampled grid removed in 2018. `DP1` now matches the released 3.6.x MEX, Python `DP` and R `DPo` to machine precision

# 3.7.0
* add warping penalties (`none`, `roughness`, `l2gam`, `l2psi`, `geodesic`) to `DynamicProgrammingQ2`
* `optimum_reparam` takes a `penalty` argument, honored by the `DP`, `DP1` and `RBFGS` methods; when omitted each method keeps its previous penalty (`none` for `DP1`, `roughness` for `DP` and `RBFGS`)
* replace vendored minFunc with Optimization Toolbox `fminunc`
* replace vendored `LMFnlsq` sphere fitting with Optimization Toolbox `lsqnonlin`
* replace vendored fdaM `basis` package with Curve Fitting Toolbox `spcol`/`augknt`
* replace vendored `bspline_tools` with Curve Fitting Toolbox `csapi`/`fnder`
* replace `varycolor` with the built-in `turbo` colormap
* remove dead code and strip cached outputs from `example.m`
* `buildtool` no longer builds the minFunc MEX files

# 3.6.11
* bug fixes to TB code from refactor
* update scale in curve code
* bugfixes

# 3.6.10
* bugfixes

# fdasrvf 3.6.9
* update armadillo
* update plotting colors to match other packages
* add project to fpca functions
* add h transform to fpca functions for warping functions
* change jfpca to approach of Wu
* add var_exp to fpca functions

# fdasrvf 3.6.8
* move to buildfile, no longer run `setup.m`, now run buildtool to create mex files and test
* move to MATLAB project structure
* create toolbox build and deploy for all 3 OS to github release

# fdasrvf 3.6.7
* bugfixes, updates for MATLAB changes

# fdasrvf 3.6.6
* update package structure

# fdasrvf 3.6.5
* bugfixes and release
