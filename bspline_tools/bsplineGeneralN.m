function y=bsplineGeneralN(x,n)
% B-spline function of order n
%
% NOTE: quite inefficient implementation

y = zeros(size(x));
for k = 0:(n+1)
    y=y+(1/factorial(n))*nchoosek(n+1,k)*(-1)^k*powplus(x-k+(n+1)/2,n);
end

end

function y=powplus(x,n)

    y = zeros(size(x));
    pidx = x>=0;
    y(pidx)=x(pidx).^n;
end
    
    