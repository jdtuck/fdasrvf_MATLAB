function gam = rgam(N, sigma, num)
% RGAM Generate random warping functions
% -------------------------------------------------------------------------
% Generates random warping functions
%
% Usage: gam = rgam(N, sigma, num)
% 
% Input:
% N: length of warping function
% sigma: variance of warping functions
% num: number of warping functions
% 
% Output:
% gam: matrix of warping functions

gam = zeros(num,N);
time = linspace(0,1,N);
omega = (2*pi);
for k = 1:num
    v = zeros(1,N);
    cnt = 1;
    for l = 2:3
        alpha_i = normrnd(0,sigma);
        if mod(l,2)~=0 %odd
            v = v + alpha_i*sqrt(2).*cos(cnt.*omega.*time);
            cnt = cnt + 1;
        end
        if mod(l,2)==0 %even
            v = v + alpha_i*sqrt(2).*sin(cnt.*omega.*time);
        end
    end
    psi = exp_map(ones(1,N), v);
    gam(k,:) = cumtrapz(time,psi.^2);
end

end
