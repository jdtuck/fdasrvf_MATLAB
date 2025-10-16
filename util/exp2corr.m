function out = exp2corr(sigma2,phi,ds)
% squared-exponential correlation function
out = sigma2 * exp(-ds.^2/(2*phi^2));