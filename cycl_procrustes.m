function [Xmtch, Xp_star, p_star, O_star, obj] = cycl_procrustes(X,Y)
%CYCL_PROCRUSTES Align two matrices by a cyclic row shift and an orthogonal 
% alignment
%
%   [XMTCH, XP_STAR, P_STAR, O_STAR, OBJ] = CYCL_PROCRUSTES(X,Y)
%
%   INPUT
%       X  – n‑by‑m matrix to be shifted (rows are ordered observations).  
%       Y  – n‑by‑m reference matrix (same size as X).  
%
%   OUTPUT
%       Xmtch   – X after the optimal shift (p_star) and the optimal
%                 orthogonal rotation (O_star).  
%       Xp_star – X after the optimal shift, before rotation.  
%       p_star  – Index of the optimal cyclic shift (1 ≤ p_star ≤ n).  
%       O_star  – Orthogonal matrix (m‑by‑m) that best aligns Xp_star with Y.  
%       obj     – n‑by‑1 vector of the objective value for every possible
%                 shift; obj(p_star) is the maximal value.  
%
%   THEORETICAL BACKGROUND
%       For a given pair of matrices A and B the function
%           fobj(A,B) = sum( svd( A' * B ) )
%       is the Schatten‑1 (nuclear) norm of A'B.  Maximising this quantity
%       over orthogonal matrices O is the classic orthogonal Procrustes
%       problem, whose solution is O = U V' where
%           [U,~,V] = svd( A' * B ).
%
%       Here the rows of X may be mis‑aligned by an unknown cyclic offset.
%       The algorithm evaluates all n possible offsets, computes the
%       Schatten‑1 inner product for each, selects the offset that gives the
%       largest value, and finally solves the Procrustes problem for that
%       offset.
%
%   EXAMPLE
%       % Generate synthetic data
%       rng(1);
%       X = randn(100,3);
%       R = orth(randn(3));           % random rotation
%       Y = X * R;                    % perfectly rotated copy
%       Y = [Y(43:end,:); Y(1:42,:)]; % introduce a cyclic shift of 42 rows
%
%       % Recover shift and rotation
%       [Xmtch, Xp_star, p_star, O_star, obj] = cycl_procrustes(X,Y);
%       fprintf('Optimal permutation index: row %d (should be 42)\n', p_star);
%       disp('Rotation error (Frobenius norm):');
%       disp(norm(O_star - R,'fro'));
%
%   REFERENCES
%       * Gower, J. C., & Dijksterhuis, G. B. (2004). *Procrustes Problems*.
%       * Bhatia, R. (1997). *Matrix Analysis*.
%
%   SEE ALSO svd, orth, procrustes
%   -------------------------------------------------------------------------
%   Author : Zach Grey
%   Created: 2023‑08‑15
%   Updated: 2025‑11‑28
%   Version: 1.2

%% ---- Input validation -------------------------------------------------
if nargin < 2
    error('Two input matrices (X,Y) are required.');
end
[nX,mX] = size(X);
[nY,mY] = size(Y);
if nX ~= nY, error('X and Y must have the same number of rows.'); end
if mX ~= mY, error('X and Y must have the same number of columns.'); end
n = nX;                                 % number of possible shifts

%% ---- Objective function (Schatten‑1 norm) -----------------------------
fobj = @(A,B) sum(svd(A'*B));

% ProjY = (eye(n) - Y* ((Y'*Y) \ Y')); fobj = @(X,Y) -norm(ProjY*X,'fro'); % explicit variable projection
% fobj = @(X,Y) norm(Y'*(X*X')*Y,'fro')^2; % explicit Frobenius norm squared

% Interesting conditional equivalences: (sought Corollary)
% fobj = @(X,Y) norm(Y'*X,'fro')^2; % Shcatten lower bound (Frobenius squared)
% fobj = @(X,Y) trace(Y'*(X*X')*Y); % explicit Frobenius norm squared
% fobj = @(X,Y) max(svd(Y'*X))*sum(svd(Y'*X)); % bounded composite
% fobj = @(X,Y) sqrt(2)*max(svd(Y'*X))*norm(Y'*X,'fro'); % CS upper

%% ---- Exhaustive (brute force) search over all cyclic shifts ------------------------
VPfunc = zeros(n,1);                    % objective for each shift

parfor p = 1:n-1                        % parallel loop (parfor)
    Xp = [X(p+1:end,:); X(1:p,:)];      % cyclically shifted X
    VPfunc(p) = fobj(Xp,Y);            % evaluate objective
end
VPfunc(n) = fobj(X,Y);                  % case of no shift (p = n)

% ---- Identify optimal shift -------------------------------------------
[~,p_star] = max(VPfunc);               % index of maximal objective

% ---- Re‑construct shifted matrix --------------------------------------
Xp_star = [X(p_star+1:end,:); X(1:p_star,:)];

% ---- Solve orthogonal Procrustes for the chosen shift ----------------
[U_star,~,V_star] = svd(Xp_star'*Y);
O_star = U_star * V_star';

% ---- Apply rotation ----------------------------------------------------
Xmtch = Xp_star * O_star;

% ---- Return full objective vector --------------------------------------
obj = VPfunc;

end