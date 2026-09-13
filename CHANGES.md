# 3.6.X
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
