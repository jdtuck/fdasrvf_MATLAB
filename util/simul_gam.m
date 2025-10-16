function gam = simul_gam(u,g1,g2,t,s1,s2,tt)
% SIMUL_GAM Align two function using simultaneous alignments.
% -------------------------------------------------------------------------
%
% Usage: gam = simul_gam(u,g1,g2,t,s1,s2,tt)
%
% Input:
%   g1,g2,s1,s2     Output of 'simul_align' representing simultaneous
%                   alignment of two functions.
% 
%   u               Common parameter of g1,g2. Must be same length as
%                   g1,g2 and increasing with domain [0,1].
% 
%   t               Common parameter of s1,s2. Must be same length as 
%                   s1,s2 and increasing with domain [0,1]. 
%   
%   tt              Discretization on which to interpolate gammas. Can
%                   be any diffeo of [0,1]. Fine discretizations give best
%                   results.
%   
% Output:
%   gam             A single reparameterization that represents the same
%                   alignment as the input simultaneous alignment. Given
%                   two functions, f1,f2, composing (t,f2) with (tt,gam) 
%                   has the same effect (within some numerical error) as
%                   composing (s1,f1) with (u,g1) and composing (s2,f2)
%                   with (u,g2).
ss = tt;
gs1 = interp1_flat(u,g1,ss);
gs2 = interp1_flat(u,g2,ss);

gt1 = interp1_flat(s1,t,gs1);
gt2 = interp1_flat(s2,t,gs2);

gam = interp1_flat(gt1,gt2,tt);
end

function yy = interp1_flat(x,y,xx)
% INTERP1_FLAT Flat linear interpolation
flat = find(diff(x)==0);
n = length(flat);

if n==0
    yy = interp1(x,y,xx);
else
    
    yy = zeros(size(xx));
    
    i1 = 1;
    if flat(1)==1
        i2 = 1;
        j = xx==x(i2);
        yy(j) = min(y(i2:i2+1));
    else
        i2 = flat(1);
        j = xx>=x(i1) & xx<=x(i2);
        yy(j) = interp1(x(i1:i2),y(i1:i2),xx(j));
        i1 = i2;
    end
    for k=2:n
        i2 = flat(k);
        if i2>i1+1
            j = xx>=x(i1) & xx<=x(i2);
            yy(j) = interp1(x(i1+1:i2),y(i1+1:i2),xx(j));
        end
        j = xx==x(i2);
        yy(j) = min(y(i2:i2+1));
        i1 = i2;
    end
    i2 = length(x);
    j = xx>=x(i1) & xx<=x(i2);
    if i1+1==i2
        yy(j) = y(i2);
    else
        yy(j) = interp1(x(i1+1:i2),y(i1+1:i2),xx(j));
    end
end
end