clc;clear
time = linspace(0,1,101);
binsize = mean(diff(time));
M = length(time);
N = 40;
lam = 0.001;
center = [.3 .35];
center2 = [4 4];
sd1 = .1;
gam_sd = .3;
num_comp = 5;
f_orig = zeros(M,N*length(center));
omega = 2*pi;
cnt = 1;
for ii=1:length(center)
    tmp = normpdf(time,center(ii),.1);
    for jj = 1:N
        f_orig(:,cnt) = normrnd(center2(ii), sd1) * tmp;
        cnt = cnt + 1;
    end
end

q_orig = f_to_srvf(f_orig,time);
y_orig = ones(size(q_orig,2),1);
y_orig(N+1:2*N) = -1;

f = zeros(M, size(f_orig,2));
q = zeros(size(f));
gam_orig = rgam(M, gam_sd, 2*N).';
for jj = 1:2*N
    f(:,jj) = interp1(time, f_orig(:,jj), (time(end)-time(1)).*gam_orig(:,jj) + time(1))';
    q(:,jj) = gradient(f(:,jj), binsize)./sqrt(abs(gradient(f(:,jj), binsize))+eps);
end
option.parallel = 0;
option.closepool = 0;
option.smooth = 0;
option.sparam = 25;
option.B = [];
option.df = 20;
option.max_itr = 20;
out = elastic_logistic(f, y_orig, time, option);
