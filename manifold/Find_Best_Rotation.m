function [q2new,R] = Find_Best_Rotation(q1,q2)
% FIND_BEST_ROTATION Find best rotation between two curves
% -------------------------------------------------------------------------
% Find best rotation between two curves
% 
% Usage: [q2new,R] = Find_Best_Rotation(q1,q2)
%
% This function finds the best rotation between two curves
%
% Input:
% q1: matrix (n,T) defining T points on n dimensional curve
% q2: matrix (n,T) defining T points on n dimensional curve
% 
% Output:
% q2new: rotated new curve
% R: rotation matrix
[n,~] = size(q1);
A = q1*q2';
[U,~,V] = svd(A);
if det(A) > 0
    S = eye(n);
else
    S = eye(n);
    S(:,end) = -S(:,end);
end
R = U*S*V';
q2new = R*q2;