function [f1_align, f2_align] = pairwise_align(f1,f2,time)
addpath(genpath('DP'))
q1 = f_to_srvf(f1,time); 
q2 = f_to_srvf(f2,time);

if iscolumn(time)
    time = time';
end

if iscolumn(f1)
    f1 = f1';
    q1 = q1';
end

if iscolumn(f2);
    f2 = f2';
    q2 = q2';
end

[G,T] = DynamicProgrammingQ2(q1/norm(q1),time,q2/norm(q2),time,time,time,0);
gam0 = interp1(T,G,time);
gam = (gam0-gam0(1))/(gam0(end)-gam0(1));  % slight change on scale

f1_align = f1;
f2_align = interp1(time, f2, (time(end)-time(1)).*gam + time(1));