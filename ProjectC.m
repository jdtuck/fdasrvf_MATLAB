function qnew = ProjectC(q)
% PROJECTC Project closed curve
% -------------------------------------------------------------------------
% Project closed curve
% 
% Usage: qnew = ProjectC(q)
%
% This function projects a curve on the tangent space
%
% Input:
% q: matrix (n,T) defining T points on n dimensional curve
% 
% Output:
% qnew: matrix of SRVF
[n,T] = size(q);
if(n == 2)
    dt = 0.35;
end
if(n == 3)
    dt = 0.2;
end
epsilon = 1e-6;

iter = 1;
res = ones(1,n);
J = zeros(n,n);

s = linspace(0,1,T);

qnew = q;
qnew = qnew/sqrt(InnerProd_Q(qnew,qnew));

qnorm = zeros(1,T);
G = zeros(1,n);
C = zeros(1,300);
while (norm(res) > epsilon)
    if(iter > 300)
        break;
    end
    % Compute Jacobian
    for i = 1:n
        for j = 1:n
            J(i,j) = 3 * trapz(s,qnew(i,:) .*qnew(j,:) );
        end
    end
    J = J+ eye(size(J));
    

    for i = 1:T
        qnorm(i) = norm(qnew(:,i));
    end

    %%%%%%%%%%%%%%%%
    % Compute the residue
    for i = 1:n
        G(i) = trapz(s,qnew(i,:).*qnorm);
    end
    res = -G;

    if(norm(res) < epsilon)
        break;
    end

    x = J\(res');
    C(iter) = norm(res);
    
    delG = Basis_Normal_A(qnew);
    temp = 0;
    for i = 1:n
        temp = temp + x(i)*delG{i}*dt;
    end
    qnew = qnew + temp;
    iter = iter + 1;
end

qnew = qnew/sqrt(InnerProd_Q(qnew,qnew));
