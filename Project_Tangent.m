function fnew = Project_Tangent(f,q)
% PROJECT_TANGENT Project tangent
% -------------------------------------------------------------------------
% This function projects the tangent vector f in the
% Tangent space of {\cal C} at q
% 
% Usage:  fnew = Project_Tangent(f,q)
%
% Input:
% f: matrix (n,T) defining T points on n dimensional vector
% q: matrix (n,T) defining T points on n dimensional vector
% 
% Output:
% w_new: transmported w

[n,~] = size(q);
% Project w in T_q({\cal B}), ie the unit sphere
w = f - InnerProd_Q(f,q)*q;
% Form the basis for the Normal space of {\cal A}
g = Basis_Normal_A(q);

% Refer to the function Gram_Schmidt for the parameters
Evorth = Gram_Schmidt(g,'InnerProd_Q');
Ev = zeros(n,T,n);
% Unpack Evorth structure
for i = 1:n
    Ev(:,:,i) = Evorth{i};
end

for i = 1:n
    fnew = w - InnerProd_Q(w,Ev(:,:,i))*Ev(:,:,i);
end
