function geodmean = geodmeanSk(vdata,Error)
% geodmeanSk    geodesic mean of data on S^k (Sphere) use Log map and Exp
%               map iterated algorithm to calculate geodesic mean. 
%               Uses rotation matrix.
%
% Inputs:   vdata   - a matrix (k+1)-by-n : a column vector represents a point on S^k  
%           Error   - a scalar for optimization 
%                 
% Output:
%           geodmean    - geodesic mean on S^k
% 
%
%   See also LogNP, ExpNP, geodmeanS2q, geodmeanS2grid,rotM. 

% 2008.07.08, modified 2008.09.28
% Aug 10, 2009
% Sungkyu Jung

if nargin == 1
    Error = 1e-10;
end

vini = vdata(:,1);
% Initial candidate for geodesic mean
diff = 1;
while diff > Error
    rot= rotMat(vini);
    rotatedvdata = rot*vdata ;
    % move all data near north pole at which the candidate mean lies
    
    vnew = ExpNPd(mean(LogNPd(rotatedvdata),2));
    % project new candidate mean back to the sphere
    rotatebackvnew = rot\vnew ;
    diff = norm(rotatebackvnew - vini);
    % see the difference
    vini = rotatebackvnew;
end;

geodmean = rotatebackvnew;
