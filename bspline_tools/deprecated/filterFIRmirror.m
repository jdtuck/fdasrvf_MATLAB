function y = filterFIRmirror(b,x) 

M = (length(b)-1)/2;

xmir = mirrorpad(x, M);
y = conv(xmir,b,'valid');