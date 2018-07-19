function w_new = Parallel_Transport_C(w,q2)
% PARALLEL_TRANSPORT Parallel transport vector
% -------------------------------------------------------------------------
% This function parallel transports a vector w from
% alpha1 \in T_{\alpha1}{\cal C} to alpha2 \in T_{\alpha2}{\cal C}
% 
% Usage:  w_new = Parallel_Transport_C(w,q2)
%
% Input:
% w: matrix (n,T) defining T points on n dimensional vector
% q2: matrix (n,T) defining T points on n dimensional vector
% 
% Output:
% w_new: transmported w

lw = sqrt(InnerProd_Q(w,w));
if(lw < 0.0001)
    w_new = w;
else
    w_new = Project_Tangent(w,q2);
    w_new = w_new*lw/sqrt(InnerProd_Q(w_new,w_new));
end

return;
