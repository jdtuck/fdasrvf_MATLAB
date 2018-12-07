function rgam = randomGamma(gam,num)
% RANDOMGAMMA generative model for gamma
% -------------------------------------------------------------------------
% This function randomly generates random gamma's by imposing a probability
% model on the space of provided gammas
%
% Usage:  rgam = randomGamma(gam,num)
%
% Input:
% gam: matrix defining warping functions (NxM) - N warping functions
% num: number of samples
%
% Output:
% rgam: random warping functions

[mu,~,vec] = SqrtMean(gam);

K = cov(vec);
[U,S,~] = svd(K);
Sig = diag(S);
n = 5;
T = length(mu);
time = linspace(0,1,T);
time = time(:);

vm = mean(vec);
rgam = zeros(num, length(mu));
for k=1:num

    a = randn(1,n);
    v = zeros(size(vm));
    for i=1:n
        v = v + a(i)*sqrt(Sig(i))*U(:,i)';
    end
    psi = exp_map(mu, v);

    gam0 = cumtrapz(time,psi.^2);
    rgam(k,:) = (gam0-gam0(1))/(gam0(end)-gam0(1));  % slight change on scale

end
end