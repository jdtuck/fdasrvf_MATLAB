function [exp_inv, theta] = inv_exp_map(Psi, psi)
% INV_EXP_MAP Inverse Exponential Map
Psi = Psi(:);
psi = psi(:);

tmp = inner_product(Psi,psi);
if tmp > 1
    tmp = 1;
end
if tmp < -1
    tmp = -1;
end
theta = acos(tmp);

if (theta < 1e-10)
    exp_inv = zeros(1,length(psi));
else
    exp_inv = theta / sin(theta) * (psi - cos(theta)*Psi);
end
