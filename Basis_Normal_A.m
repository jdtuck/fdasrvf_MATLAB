function delG = Basis_Normal_A(q)
% BASIS_NORMAL_A Find Normal Basis
% -------------------------------------------------------------------------
% This function returns the vector field that forms the basis
% for normal space of {\cal A}
% 
% Usage: delG = Basis_Normal_A(q)
%
% Input:
% q: matrix (n,T) defining T points on n dimensional SRVF
% 
% Output:
% delg: Basis

[n,T] = size(q);
e = eye(n);
Ev = zeros(n,T,n);
for i = 1:n
    Ev(:,:,i) = repmat(e(:,i),1,T);
end

qnorm = zeros(1,T);
for t = 1:T
    qnorm(t) = norm(q(:,t));
end

delG = cell(1,n);
for i = 1:n
    tmp1 = repmat(q(i,:)./qnorm,n,1);
    tmp2 = repmat(qnorm,n,1);
    delG{i} = tmp1.*q + tmp2.*Ev(:,:,i);    
end