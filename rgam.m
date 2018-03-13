function gam = rgam(N, sigma, num)
% Random Warping
%
% Generates random warping functions
%
% @param N length of warping function
% @param sigma variance of warping functions
% @param num number of warping functions
% @return gam warping functions

gam = zeros(num,N);
TT = N - 1;
time = linspace(0,1,TT);
mu = sqrt(ones(1,N-1)*TT/(N-1));
omega = (2*pi)/TT;
for k = 1:num
    alpha_i = normrnd(0,sigma);
    v = alpha_i * ones(1,TT);
    cnt = 1;
    for l = 2:10
        alpha_i = normrnd(0,sigma);
        if mod(l,2)~=0 %odd
            v = v + alpha_i*sqrt(2).*cos(cnt.*omega.*time);
            cnt = cnt + 1;
        end
        if mod(l,2)==0 %even
            v = v + alpha_i*sqrt(2).*sin(cnt.*omega.*time);
        end
    end
    v = v - ((mu*v.')*mu)/TT;
    vn = norm(v)/sqrt(TT);
    psi = cos(vn).*mu + sin(vn).*v./vn;
    gam(k,:) = [0 cumsum(psi.*psi)]./N;
end

end
