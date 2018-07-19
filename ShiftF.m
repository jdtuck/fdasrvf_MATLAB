function pn = ShiftF(p,tau)
% SHIFTF Shifts curve starting point by tau
% -------------------------------------------------------------------------
% This function shifts a curves starting point by tau
%
% Usage: pn = ShiftF(p,tau)
%
% Input:
% p: matrix (n,T) defining T points on n dimensional curve
% tau: integer to shift curve
%
% Output
% pn: shifted curve
[~,T] = size(p);
if(tau == 0)
    pn = p;
    return;
end

if tau > 0
    pn(:,1:T-tau) = p(:,tau+1:T);
    pn(:,T-tau+1:T) = p(:,1:tau);
else
    t = abs(tau)+1;
    pn(:,1:T-t+1) = p(:,t:T);
    pn(:,T-t+2:T) = p(:,1:t-1);
end