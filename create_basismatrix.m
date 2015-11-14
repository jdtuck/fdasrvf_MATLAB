function basismat = create_basismatrix(x, df, degree)
addpath(genpath('basis'))
tmp = create_bspline_basis([min(x) max(x)],df,degree);
rangex   = getbasisrange(tmp);
x        = linspace(rangex(1),rangex(2),length(x))';
basismat = full(eval_basis(x, tmp));