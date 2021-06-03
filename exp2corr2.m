function out = exp2corr2(phi,ds)
% squared-exponential correlation function
out = exp(-ds.^2/(2*phi^2));
