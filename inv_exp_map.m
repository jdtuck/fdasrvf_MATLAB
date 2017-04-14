function [exp_inv, theta] = inv_exp_map(Psi, psi)

theta = acos(inner_product(Psi,psi));

if (theta < 1e-10)
    exp_inv = zeros(1,length(psi));
else
    exp_inv = theta / sin(theta) * (psi - cos(theta)*Psi);
end
