function colmap = rainbowcol(n)
%RAINBOWCOL gives rainbow color matrix of for length n,
% modified from Steve Marron's Matlab programs.

if n ==2
    colmap = [0 0 1; 1 0 0];
elseif n==3
    colmap = [0 0 1; 0.3 0.7 0.3 ; 1 0 0];
else
    colmap = colormap(jet(n));
end

    %
% nfifth = ceil((n - 1) / 5) ;
% del = 1 / nfifth ;
% vwt = (0:del:1)' ;
% colmap = [flipud(vwt), zeros(nfifth+1,1), ones(nfifth+1,1)] ;
% colmap = colmap(1:size(colmap,1)-1,:) ;
% %  cutoff last row to avoid having it twice
% colmap = [colmap; ...
%     [zeros(nfifth+1,1), vwt, ones(nfifth+1,1)]] ;
% colmap = colmap(1:size(colmap,1)-1,:) ;
% %  cutoff last row to avoid having it twice
% colmap = [colmap; ...
%     [zeros(nfifth+1,1), ones(nfifth+1,1), flipud(vwt)]] ;
% colmap = colmap(1:size(colmap,1)-1,:) ;
% %  cutoff last row to avoid having it twice
% colmap = [colmap; ...
%     [vwt, ones(nfifth+1,1), zeros(nfifth+1,1)]] ;
% colmap = colmap(1:size(colmap,1)-1,:) ;
% %  cutoff last row to avoid having it twice
% colmap = [colmap; ...
%     [ones(nfifth+1,1)], flipud(vwt), zeros(nfifth+1,1)] ;
% 
% %  note: put this together upside down
