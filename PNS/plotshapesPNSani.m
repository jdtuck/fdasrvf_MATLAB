function plotshapesPNSani(PCnum,steps,savestr,EuclidData,PNS,joinline)
% plotshapesPNSani : animated motion of shapes for PNSshape
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
if nargin < 6
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
movieframe = [m:length(steps) length(steps):-1:1 1:m];

% initialize movie frame
clf;
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
X = H'*reshape(PCvec(:,m),(k-1),2);
plot(X(joinline,1),X(joinline,2),'-r','Linewidth',3);
if ~isempty(setdiff(1:k,joinline))
    scatter(X(setdiff(1:k,joinline),1),X(setdiff(1:k,joinline),2),'or');
end
text(0,0,'1','Color','k');
hold off;axis off;
limtmp=axisSM(Xmat(:,1),Xmat(:,2));
limfix = [min(limtmp([1 3])), 1.1*max(limtmp([2 4]))];
set(gca,'Xlim',limfix,'Ylim',limfix);
f = getframe;
[im,map] = rgb2ind(f.cdata,256,'nodither');
im(1,1,1,length(movieframe)) = 0;

% capture movie frame
for i=1:length(movieframe)
    clf;hold on;
    for j=1:n;
        X = H'*reshape(preshapeT(:,j),(k-1),d);
        plot(X(joinline,1),X(joinline,2),':','Color',0.7*[1 1 1],'Marker','.');
        if ~isempty(setdiff(1:k,joinline))
            scatter(X(setdiff(1:k,joinline),1),X(setdiff(1:k,joinline),2),'MarkerEdgeColor',0.7*[1 1 1],'Marker','.');
        end
    end
    X = H'*reshape(PCvec(:,movieframe(i)),(k-1),2);
    plot(X(joinline,1),X(joinline,2),'-r','Linewidth',3);
    if ~isempty(setdiff(1:k,joinline))
        scatter(X(setdiff(1:k,joinline),1),X(setdiff(1:k,joinline),2),'or');
    end
    if strcmp(PNS.itype,'great')
        texttoprint = ['Along PNG' num2str(PCnum) ' at ' num2str(steps(movieframe(i)),'%2.2f') ' \sigma'];
    else
        texttoprint = ['Along PNS' num2str(PCnum) ' at ' num2str(steps(movieframe(i)),'%2.2f') ' \sigma'];
    end
    text(limfix(1)*0.9,limfix(2)*0.9,texttoprint,'FontSize',18,'Color','k')
    hold off;axis off;set(gca,'Xlim',limfix,'Ylim',limfix);
    %set(gca,'nextplot','replacechildren','visible','off')
    f = getframe;
    im(:,:,1,i) = rgb2ind(f.cdata,map,'nodither');
end
if strcmp(PNS.itype,'great')
    imwrite(im,map,[savestr '_PNG' num2str(PCnum) '.gif'],'DelayTime',0,'LoopCount',100) %g443800
else
    imwrite(im,map,[savestr '_PNS' num2str(PCnum) '.gif'],'DelayTime',0,'LoopCount',100) %g443800
end