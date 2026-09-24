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

[n,T] = size(q);
% Form the basis for the Normal space of {\cal A} and orthonormalize it
% together with q so that the projection removes both the radial
% direction of the unit sphere {\cal B} and the closure normals
g = Basis_Normal_A(q);

% Refer to the function Gram_Schmidt for the parameters
Evorth = Gram_Schmidt([{q}, g],'InnerProd_Q');
nb = length(Evorth);
Ev = zeros(n,T,nb);
% Unpack Evorth structure
for i = 1:nb
    Ev(:,:,i) = Evorth{i};
end

fnew = f;
for i = 1:nb
    fnew = fnew - InnerProd_Q(f,Ev(:,:,i))*Ev(:,:,i);
end
