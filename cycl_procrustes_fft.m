function [Xmtch, Xp_star, p_star, O_star, obj] = cycl_procrustes_fft(X,Y)
%CYCL_PROCRUSTES_FFT  Fast cyclic Procrustes using FFT
%
%   [Xmtch, Xp_star, p_star, O_star, obj] = cycl_procrustes_fft(X,Y)
%
%   INPUT
%       X : n‑by‑m matrix whose rows may be cyclically shifted
%       Y : n‑by‑m reference matrix (same size as X)
%
%   OUTPUT
%       Xmtch   : X after optimal shift AND optimal orthogonal rotation
%       Xp_star : X after optimal shift, *before* rotation
%       p_star  : index (0‑based, i.e. number of rows moved from the front)
%                 of the optimal cyclic shift
%       O_star  : orthogonal m‑by‑m matrix that aligns Xp_star with Y
%       obj     : n‑by‑1 vector containing the nuclear‑norm objective for
%                 every possible shift; obj(p_star+1) is the maximal value.
%
%   NOTE
%       The implementation works for any m≥1, but it shines when m is
%       small (2–3) and n is large (10⁴–10⁶).  The algorithm is
%       O(m² n log n) because the cross‑correlations are evaluated with a
%       single FFT per column‑pair instead of a full O(n) matrix product
%       for each shift.
%
%   Author : Zach Grey (adapted by OpenAI ChatGPT)
%   Date   : 2025‑12‑01
%   ---------------------------------------------------------------

% -----------------------------------------------------------------
% 0. Input validation (mirrors the original routine)
% -----------------------------------------------------------------
if nargin < 2
    error('Two input matrices (X,Y) are required.');
end
[ nX , mX ] = size(X);
[ nY , mY ] = size(Y);
if nX ~= nY, error('X and Y must have the same number of rows.'); end
if mX ~= mY, error('X and Y must have the same number of columns.'); end
n = nX;                 % number of possible cyclic shifts
m = mX;

% -----------------------------------------------------------------
% 1. Compute all circular cross‑correlations Xj ↔ Yk  (FFT)
% -----------------------------------------------------------------
% The FFT length is exactly n because we want a *circular* correlation.
% (If n is not a power of two, MATLAB's FFT is still O(n log n).)
Xf = fft(X, n, 1);      % n‑by‑m, FFT along rows
Yf = fft(Y, n, 1);      % n‑by‑m

% Allocate a 3‑D array that will hold C(p) = X(p)ᵀY for all p.
% The layout will be (p , j , k) → C(p,j,k) .
C = zeros(n, m, m);     % n shifts, each an m×m matrix

% Fill C by looping over column pairs (j,k).  Because m is tiny this loop
% is negligible compared with the FFTs.
for j = 1:m
    % Element‑wise product of the j‑th column of X with *all* columns of Y
    prod_jk = Xf(:,j) .* conj(Yf);            % n‑by‑m (complex)
    % Inverse FFT gives the circular cross‑correlation for every shift.
    % Real part is taken (imaginary part is ≈ 0 due to round‑off).
    corr_jk = real(ifft(prod_jk, [], 1));      % n‑by‑m
    % Store: C(p,j,k) = corr_jk(p,k)
    C(:,j,:) = corr_jk;                       % reshape automatically
end

% At this point C(p,:,:) == X(p)ᵀY for every cyclic shift p (0‑based).

% -----------------------------------------------------------------
% 2. Evaluate the nuclear‑norm objective for each shift
% -----------------------------------------------------------------
obj = zeros(n,1);
for p = 1:n
    % Extract the m×m matrix for shift p‑1 (MATLAB indices are 1‑based)
    Cp = squeeze(C(p,:,:));      % m‑by‑m
    % Nuclear norm = sum of singular values (Schatten 1-norm)
    s  = svd(Cp);
    obj(p) = sum(s);
end

% -----------------------------------------------------------------
% 3. Locate the optimal shift
% -----------------------------------------------------------------
[~,p_star] = max(obj);           % p_star is 1‑based (MATLAB‑style)

% -----------------------------------------------------------------
% 4. Build X shifted by the optimal amount (no rotation yet)
% -----------------------------------------------------------------
Xp_star = [ X(p_star:end,:) ; X(1:p_star-1,:) ];   % size n×m

% -----------------------------------------------------------------
% 5. Orthogonal Procrustes for the chosen shift
% -----------------------------------------------------------------
[U,~,V] = svd( Xp_star.' * Y );
O_star  = U * V.';               % orthogonal m‑by‑m matrix

% -----------------------------------------------------------------
% 6. Apply the rotation
% -----------------------------------------------------------------
Xmtch = Xp_star * O_star;

% -----------------------------------------------------------------
% 7. Return (MATLAB‑style 1‑based index for p_star)
% -----------------------------------------------------------------
% The output ``obj`` already contains the full objective vector.
end