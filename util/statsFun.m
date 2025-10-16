function out = statsFun(mat)
out = quantile(mat.',[0.025,.975]);
end