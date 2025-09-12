function xy = stereographicproj(data,sp);
[d n] = size(data);

if d ~= 3
    return;
end

rotdata = rotMat(sp,[0,0,-1]')*data;
xy = [rotdata(1,:)./(1-rotdata(3,:));
      rotdata(2,:)./(1-rotdata(3,:))];