% simul_gam -- align one function to another given two simultaneous
%           alignments.
%
% Arguments:
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
% Return value:
%   gam             A single reparameterization that represents the same
%                   alignment as the input simultaneous alignment. Given
%                   two functions, f1,f2, composing (t,f2) with (tt,gam) 
%                   has the same effect (within some numerical error) as
%                   composing (s1,f1) with (u,g1) and composing (s2,f2)
%                   with (u,g2).
function gam = simul_gam(u,g1,g2,t,s1,s2,tt)
    
    ss = tt;
    gs1 = interp1(u,g1,ss);
    gs2 = interp1(u,g2,ss);
    
    gt1 = interp1_flat(s1,t,gs1);
    gt2 = interp1_flat(s2,t,gs2);
    
    gam = interp1_flat(gt1,gt2,tt);
end