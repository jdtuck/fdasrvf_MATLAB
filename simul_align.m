% simul_align -- align two functions by simultaneous reparameterization
% 
% Arguments:
%   f1,f2       Data functions to be aligned. They should be the same 
%               length but no special assumption is needed about their 
%               respective parameterizations. It's only assumed that 
%               they are strictly increasing with domain [0,1].
% Return values:
%   s1,s2       Arc-length parameterizations of f1,f2. I.e., the pairs
%               of vectors (s1,f1) and (s2,f2) represent sample points 
%               of arc-length parameterized functions on [0,1].
%
%   g1,g2       Reparameterizations to simultaneously align f1 and f2.
%               The two g functions are represented by vectors of the same
%               length and are assumed to have a common parameter. The
%               common parameter isn't given here but it need only be a
%               vector of the same length which is increasing and has the
%               domain [0,1].
%
%   ext1,ext2   Indices of local extrema in f1,f2
%
%   mpath       Indices of matched extrema in ext1,ext2

function [s1,s2, g1,g2, ext1,ext2, mpath] = simul_align(f1,f2)
    
    % parameterize by arc-length
    s1 = arclength(f1);
    s2 = arclength(f2);
    
    len1 = max(s1);
    len2 = max(s2);
    
    f1 = f1/len1;
    s1 = s1/len1;
    f2 = f2/len2;
    s2 = s2/len2;
    
    % get srvf (should be +/-1)
    q1 = diff(f1)./diff(s1);
    q1(diff(s1)==0) = 0;
    q2 = diff(f2)./diff(s2);
    q2(diff(s2)==0) = 0;
    
    % get extreme points
    [ext1, d1] = extrema_1s(s1,q1);
    [ext2, d2] = extrema_1s(s2,q2);
%         sprintf('segment %i',i)
    [~,~,mpath] = match_ext(s1,ext1,d1,s2,ext2,d2);
    
    te1 = s1(ext1);
    te2 = s2(ext2);
        
    [g1,g2] = simul_reparam(te1,te2,mpath);

end