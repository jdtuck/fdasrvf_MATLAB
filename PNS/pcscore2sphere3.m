function S_star = pcscore2sphere3(n_pc, X_hat, Xs, Tan, V)
    %  Converts principal component scores to points on the sphere.

    % Usage: pcscore2sphere3(n_pc, X_hat, Xs, Tan, V)

    % n_pc: number of principal components
    % X_hat: (d x n) matrix of principal component scores
    % Xs: (d x n) matrix of original data
    % Tan: tangent space at the origin
    % V: rotation matrix to align the tangent space with the sphere
    % return: (d x n) matrix of points on the sphere

    (d, n) = size(Tan);
    W = zeros(d,n);
    for i = 1:n
        W(:, i) = (acos(sum(Xs(i,:) .* Xhat)) .* Tan(:,i) / norm(Tan(:,i)));
    end
    
    lam = zeros(n,d);
    for i = 1:n
        for j = 1:n_pc
            lam(i,j) = sum(W(:,i) .* V(:,j));
        end
    end

    U = zeros(n,d);
    for i = 1:n
        for j = 1:n_pc
            U(i,:) = U(i,:) + lam(i,j) * V(:,j);
        end
    end
    S_star = zeros(n, n_pc+1);
    for i=1:n
        U_norm = norm(U(i,:));
        S_star(i, :) = [sin(U_norm) / U_norm * lam(i, 1:n_pc), 0, cos(U_norm)];
    end