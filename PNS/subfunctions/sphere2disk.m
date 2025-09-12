function xy = sphere2disk(data,sp);
% Lambert azimuthal equal-area projection
[d n] = size(data);

if d ~= 3
    return;
end

rotdata = rotMat(sp,[0,0,-1]')*data;

xy = [rotdata(1,:).*sqrt(2./(1-rotdata(3,:)));
      rotdata(2,:).*sqrt(2./(1-rotdata(3,:)))];
  
% equators
% theta = linspace(0,2*pi,40);
% plot(cos(theta)*sqrt(1), sin(theta)*sqrt(1),':k')