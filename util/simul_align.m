function [s1,s2, g1,g2, ext1,ext2, mpath] = simul_align(f1,f2)
% SIMUL_ALIGN Align two functions by simultaneous reparameterization
% -------------------------------------------------------------------------
%
% Usage: [s1,s2, g1,g2, ext1,ext2, mpath] = simul_align(f1,f2)
%
% Input:
%   f1,f2       Data functions to be aligned. They should be the same
%               length but no special assumption is needed about their
%               respective parameterizations. It's only assumed that
%               they are strictly increasing with domain [0,1].
%
% Output:
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

function t1 = arclength(f)
% ARCLENGTH Calculate arclength
t1 = zeros(length(f),1);
t1(1) = 0;
t1(2:end) = abs(diff(f));
t1 = cumsum(t1);
end

function [ext, d] = extrema_1s(t,q)
% EXTREMA_1s Finds extrema in SIMUL Warping
q = round(q);

if q(1)~=0
    d = -q(1);
else
    d = q(q~=0);
    d = d(1);
end

ext = find(diff(q))+1;

ext(2:end+1) = ext;
ext(1) = 1;
ext(end+1) = length(t);
end


function [g1, g2] = simul_reparam(te1, te2, mpath)
g1(1) = 0;
g2(1) = 0;

if mpath(1,2) == 2
    g1 = [g1; 0];
    g2 = [g2; te2(2)];
elseif mpath(1,1) == 2
    g1 = [g1; te1(2)];
    g2 = [g2; 0];
end

m = size(mpath,1);
for i=1:m-1
    [gg1,gg2] = simul_reparam_segment(mpath(i,:), mpath(i+1,:), te1, te2);
    
    g1 = [g1; gg1];
    g2 = [g2; gg2];
end

n1 = length(te1);
n2 = length(te2);
if (mpath(end,1) == n1-1) || (mpath(end,2) == n2-1)
    g1 = [g1; 1];
    g2 = [g2; 1];
end
end

function [gg1,gg2] = simul_reparam_segment(src, tgt, te1, te2)
i1 = src(1)+1:2:tgt(1);
i2 = src(2)+1:2:tgt(2);

v1 = sum(te1(i1)-te1(i1-1));
v2 = sum(te2(i2)-te2(i2-1));
R = v2/v1;

a1=src(1); a2=src(2);
t1=te1(a1); t2=te2(a2);
u1=0; u2=0;

gg1=[]; gg2=[];

while (a1<tgt(1)) && (a2<tgt(2))
    if a1==tgt(1)-1 && a2==tgt(2)-1
        a1=tgt(1); a2=tgt(2);
        gg1 = [gg1; te1(a1)];
        gg2 = [gg2; te2(a2)];
    else
        p1 = (u1 + te1(a1+1) - t1)/v1;
        p2 = (u2 + te2(a2+1) - t2)/v2;
        if p1<p2
            lam = t2 + R*(te1(a1+1)-t1);
            gg1 = [gg1; te1(a1+1); te1(a1+2)];
            gg2 = [gg2; lam; lam];
            u1 = u1 + te1(a1+1) - t1;
            u2 = u2 + lam - t2;
            t1 = te1(a1+2);
            t2 = lam;
            a1 = a1+2;
        else
            lam = t1 + (1/R)*(te2(a2+1)-t2);
            gg1 = [gg1; lam; lam];
            gg2 = [gg2; te2(a2+1); te2(a2+2)];
            u1 = u1 + lam - t1;
            u2 = u2 + te2(a2+1) - t2;
            t1 = lam;
            t2 = te2(a2+2);
            a2 = a2+2;
        end
    end
end
end

function [D,P, mpath] = match_ext(t1,ext1,d1,t2,ext2,d2)
% MATCH_EXT SIMUL matching of extents
te1 = t1(ext1);
te2 = t2(ext2);

% We'll pad each sequence to start on a 'peak' and end on a 'valley'
pad1 = [0,0];
pad2 = [0,0];
if d1==-1
    te1(2:end+1) = te1;
    te1(1) = te1(2);
    pad1(1) = 1;
end
if mod(length(te1),2)==1
    te1(end+1) = te1(end);
    pad1(2) = 1;
end

if d2==-1
    te2(2:end+1) = te2;
    te2(1) = te2(2);
    pad2(1) = 1;
end
if mod(length(te2),2)==1
    te2(end+1) = te2(end);
    pad2(2) = 1;
end

n1 = length(te1);
n2 = length(te2);

% initialize weight and path matrices
D = zeros(n1,n2);
P = zeros(n1,n2,2);

for i=1:n1
    for j=1:n2
        if mod(i+j,2)==0
            for ib = (i-1):-2:(1+mod(i,2))
                for jb = (j-1):-2:(1+mod(j,2))
                    icurr = ib+1:2:i;
                    jcurr = jb+1:2:j;
                    W = sqrt(sum(te1(icurr)-te1(icurr-1)))...
                        *sqrt(sum(te2(jcurr)-te2(jcurr-1)));
                    if ~isreal(W)
                        keyboard
                    end
                    Dcurr = D(ib,jb) + W;
                    if Dcurr > D(i,j)
                        D(i,j) = Dcurr;
                        P(i,j,:) = [ib,jb];
                    end
                end
            end
        end
    end
end

D = D(1+pad1(1):end-pad1(2), 1+pad2(1):end-pad2(2));
P = P(1+pad1(1):end-pad1(2), 1+pad2(1):end-pad2(2), :);
P(:,:,1) = P(:,:,1)-pad1(1);
P(:,:,2) = P(:,:,2)-pad2(1);

% retrieve best path

if pad1(2)==pad2(2)
    mpath = size(D);
elseif D(end-1,end) > D(end,end-1)
    mpath = size(D) - [1,0];
else
    mpath = size(D) - [0,1];
end
prev_vert = reshape(P(mpath(1,1),mpath(1,2),:),1,2);
while all(prev_vert>0)
    mpath = [prev_vert; mpath];
    prev_vert = reshape(P(mpath(1,1),mpath(1,2),:),1,2);
end

end

