function EuclidData = fastPNSs2e(spheredata, PNS)
% Fast PNS coordinate transform from Sphere to Euclidean-type residual \

% Usage: EuclidData = fastPNSe2s(spheredata, PNS)

%  where 'res' is d x m data matrix in PNS coordinate system (for any
%  m >= 1), 'PNS' is the structural array

muhat = PNS.muhat;
pca = PNS.pca;
n_pc = PNS.n_pc;

Xs = spheredata';
n = size(Xs,1);
for i = 1:n
    Xs(i,:) = Xs(i, :) / norm(Xs(i, :));
end

TT = Xs;
for i = 1:n
    TT(i,:) = Xs(i,:) - sum(Xs(i,:) .* muhat) .* muhat;
end

spdata = pcscore2sphere3(n_pc, muhat, Xs, TT', pca);

EuclidData = PNSs2e(spdata', PNS);
