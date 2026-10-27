function [exp_inv, theta] = inv_exp_map(Psi, psi)
% INV_EXP_MAP Inverse Exponential Map
% The angle is computed with atan2 of the tangent and normal components
% rather than acos of the inner product, which loses about half the
% significant digits for nearby points
Psi = Psi(:);
psi = psi(:);

ip = inner_product(Psi,psi);
w = psi - ip*Psi;
wnorm = L2norm(w);
theta = atan2(wnorm, ip);

if (theta < 1e-10)
    exp_inv = zeros(length(psi),1);
else
    exp_inv = theta / wnorm * w;
end
