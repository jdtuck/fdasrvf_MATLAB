function plotshapesPNS(PCnum,steps,EuclidData,PNS,joinline)
% plotshapesPNS : Draw Principal Arc trajectory on current axis
%
% Inputs:
% PCnum : chose Principal component number (default = 1)
% steps : vector of steps along the PC in the unit of st.dev
%        (default = -2:0.25:2)
% savestr : string for .gif file name
% joinline: a vector of numbers to join lines (default = 1:k)
% others are output of PNSshapes
%
%
% % example:
% steps = -2:0.25:2;
% PCnum = 1;
% savestr = 'shape';
% aniPNSshape(PCnum,steps,savestr,EuclidData,PNS)
%
% See also.
%
% Sungkyu Jung


[kk n]=size(EuclidData); k = (kk+1)/2 + 1; d = 2;
if nargin < 5
    joinline = [1:k 1];
end

% vector of st.dev
stdevPNS = sqrt(sum(abs(EuclidData).^2, 2) / n);
% matrix of direction vectors
udir = eye(2*(k-1) - 1);

H = Helmertsub(k);
ptEval =  udir(:,PCnum)*stdevPNS(PCnum)*steps ;
% evaluation points on pre-shape space
PCvec= PNSe2s(ptEval,PNS);
% samples on pre-shape space
preshapeT= PNSe2s(EuclidData,PNS);

m = (length(steps)+1)/2;
%movieframe = [m:length(steps) length(steps):-1:1 1:m];

% initialize movie frame
%clf;
hold on;
Xmat = zeros(k*n,2);
for j=1:n;
    X = H'*reshape(preshapeT(:,j),(k-1),d);
    Xmat((1:k)+(j-1)*n,:) = X;
    plot(X(joinline,1),X(joinline,2),':','Color',0.7*[1 1 1],'Marker','.');
    if ~isempty(setdiff(1:k,joinline))
        scatter(X(setdiff(1:k,joinline),1),X(setdiff(1:k,joinline),2),'MarkerEdgeColor',0.7*[1 1 1],'Marker','.');
    end
end

TrajectoryX = zeros(k,d,length(steps));
TrajDrawX = zeros(length(steps),d,k);
for j = 1:length(steps) 
    X = H'*reshape(PCvec(:,j),(k-1),2);
    TrajectoryX(:,:,j) = X;
    for i = 1:k
    TrajDrawX(j,:,i) = X(i,:);
    end
end
    
scatter(TrajectoryX(:,1,m),TrajectoryX(:,2,m),...
    'o','filled','MarkerEdgeColor',[0.8 0 0.8],'MarkerFaceColor',[0.8 0 0.8]); % PNSmean
for i = 1:k
   plot(TrajDrawX(1:m,1,i),TrajDrawX(1:m,2,i),'-r','Linewidth',1.5);
   plot(TrajDrawX(m:end,1,i),TrajDrawX(m:end,2,i),'-b','Linewidth',1.5);
end

hold off;

limtmp=axisSM(Xmat(:,1),Xmat(:,2));
axis equal
set(gca,'Xlim',limtmp(1:2),'Ylim',limtmp(3:4));
% 
% limtmp=axisSM(Xmat(:,1),Xmat(:,2));
% limfix = [min(limtmp([1 3])), 1.1*max(limtmp([2 4]))];
% set(gca,'Xlim',limfix,'Ylim',limfix);