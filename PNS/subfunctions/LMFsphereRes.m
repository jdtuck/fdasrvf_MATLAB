function rr = LMFsphereRes(center)
% LMFSPHERERES : Auxiliary function for LMFsphereFit.
%              -Calculates residuals of circle fit for the given center.
%
%   See also LMFsphereFit, LMFnlsq.

% Last updated Aug 10, 2009
% Sungkyu Jung

global data
global greatCircle
di = sqrt(sum(...
    (data - repmat(center,1,size(data,2))).^2 ...
    ));
r = pi/2;
if greatCircle == 0;
    r= mean(di);
end
rr = (di-r)';

