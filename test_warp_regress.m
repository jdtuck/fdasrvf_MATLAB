clc;clear
time = linspace(0,1,101);
binsize = mean(diff(time));
M = length(time);
N = 30;
lam = 0.001;
center = [.35 .5 .65];
center2 = [4 3.7 4];
sd1 = .05;
gam_sd = .25;
num_comp = 5;
f_orig = zeros(M,N*length(center));
omega = 2*pi;
cnt = 1;
for ii=1:length(center)
    tmp = normpdf(time,center(ii),.075);
    for jj = 1:N
        f_orig(:,cnt) = normrnd(center2(ii), sd1) * tmp;
        cnt = cnt + 1;
    end
end

q_orig = f_to_srvf(f_orig,time);
y_orig = zeros(size(q_orig,2),1);
alpha_t = 0;
b1 = sin(omega*time);
b2 = cos(omega*time);
B = [b1' b2'];
bt= .5*b1 + .9*b2;
cnt = 1;
for ii = 1:length(center)
    tmp = normpdf(time,center(ii),0.075);
    tmp2 = f_to_srvf((normrnd(center2(ii),sd1)*tmp)',time).';
    for jj = 1:N
        y_orig(cnt) = alpha_t + trapz(time, tmp2.*bt) + normrnd(0, 0.01);
        cnt = cnt + 1;
    end
end

f = zeros(M, size(f_orig,2));
q = zeros(size(f));
cnt = 1;
for ii = 1:length(center)
    gam_orig = rgam(M, gam_sd, N).';
    for jj = 1:N
        f(:,cnt) = interp1(time, f_orig(:,cnt), (time(end)-time(1)).*gam_orig(:,jj) + time(1))';
        q(:,cnt) = gradient(f(:,cnt), binsize)./sqrt(abs(gradient(f(:,cnt), binsize))+eps);
        cnt = cnt + 1;
    end
end
option.parallel = 0;
option.closepool = 1;
option.smooth = 0;
option.sparam = 25;
option.B = B;
option.df = 20;
option.max_itr = 20;
out = elastic_regression(f, y_orig, time,0, option);
