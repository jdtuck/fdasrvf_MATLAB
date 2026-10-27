function [exp_inv, theta] = inv_exp_map(Psi, psi)
% INV_EXP_MAP Inverse Exponential Map
Psi = Psi(:);
psi = psi(:);

% angle from chord lengths: for unit vectors ||psi-Psi|| = 2sin(theta/2) and
% ||psi+Psi|| = 2cos(theta/2). Unlike acos(<Psi,psi>), this keeps full
% relative precision for nearby points (acos floors theta at ~1e-8)
theta = 2*atan2(L2norm(psi-Psi), L2norm(psi+Psi));

if (theta < 1e-10)
    exp_inv = zeros(1,length(psi));
else
    exp_inv = theta / sin(theta) * (psi - cos(theta)*Psi);
end
