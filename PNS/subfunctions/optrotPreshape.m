function [rotz M]=optrotPreshape(z,mu)
% optrotPreshape optimal rotation of preshape z onto mu
% returns the rotated z1 (rotz1) and the optimal rotation 
% matrix M such that if z1, z2 are (k-1) x 2 preshape matrix,
% then M minimizes || mu - z M ||.

[vv dd uu] = svd(mu'*z);
refelctionBestflag = det(mu'*z)<0;
if refelctionBestflag
    % then reflection gives the optimal fit
    % thus force it to be a rotation
    uu(:,end) = -uu(:,end);
    %dd(end,end) = -dd(end,end);
    %disp(num2str(det(mu'*z)));
end

M = uu*vv';
rotz = z*M;

criteria = mu'*rotz - rotz'*mu;
if abs(criteria(2,2)) > 1e-10;
    disp('warning: rotation is not valid')
end