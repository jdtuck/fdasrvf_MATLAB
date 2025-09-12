function plotshapes(varargin)
% PLOTSHAPES shape configuration overlay
% plotshapes(A), with k x 2 x n array of configuration matrix A
% plotshapes(A,joinline), joinline: a vector of indices to connect (default [1:k 1])
% plotshapes(A,joinline,paramstruct), see below for paramstruct.
% plotshapes(resmat,PNS), with resmat, PNS are output from PNSshape.m
% plotshapes(resmat,PNS,joinline)
% plotshapes(resmat,PNS,joinline,paramstruct)
% 
% Input: paramstruct - a Matlab structure of input parameters
%                    Use: "help struct" and "help datatypes" to
%                         learn about these.
%                    Create one, using commands of the form:
%
%     paramstruct = struct('field1',values1,...
%                          'field2',values2,...
%                          'field3',values3) ;
%
%                          where any of the following can be used,
%                          these are optional, misspecified values
%                          revert to defaults
%
%    fields           values
%
%    titlestring      string for title
%
%
%    icolor           0  fully black and white version (everywhere)
%                     string (any of 'r', 'g', 'b', etc.) that single color
%                     1  (default)  color version (Matlab 7 color default)
%                     2  time series version (ordered spectrum of colors)
%                     nx3 color matrix:  a color label for each data point
%                             (to be used everywhere, except SiZer & QQ
%                              useful for comparing classes)
%
%    markerstr        Can be either a single string with symbol to use for marker,
%                         e.g. 'o' (default), '.', '+', 'x'
%                         (see "help plot" for a full list)
%                     Or a character array (n x 1), of these symbols,
%                         One for each data vector, created using:  strvcat
%


if nargin==0, help plotshapes, return, end     %   Display help
data = varargin{1};
isConfig = length(size(data))==3 || ~iscell(varargin(2));

switch isConfig
    case true % then input is array of configuration matrices
        [k d n]=size(data);
        
        if nargin < 2
            joinline = [1:k 1];
        else
            joinline = varargin{2};
        end
        
        icolor = repmat([0 0 0],n,1);
        markerstr = repmat('.',n,1);
        titlestring = 'shape configuration';
        
        if nargin > 2
            paramstruct = varargin{3};
            if isfield(paramstruct,'titlestring') ;    %  then change to input value
                titlestring = paramstruct.titlestring;
            end ;
            if isfield(paramstruct,'icolor') ;    %  then change to input value
                icolor = paramstruct.icolor;
                if length(icolor)==1 && icolor <  4
                    icolor = repmat('k',n,1);
                end
                if ischar(icolor) || size(icolor,1) == 1
                    icolor = repmat(icolor,n,1);
                end
                
            end ;
            if isfield(paramstruct,'markerstr') ;    %  then change to input value
                markerstr = paramstruct.markerstr;
                if ischar(markerstr)
                    markerstr = repmat(markerstr,n,1);
                end
            end ;
        end

        hold on;
        Xmat = zeros(k*n,2);
        for i=1:n;
            X = data(:,:,i);
            Xmat((1:k)+(i-1)*n,:) = X;
            plot(X(joinline,1),X(joinline,2),':','Color',icolor(i,:),'Marker',markerstr(i))
            if ~isempty(setdiff(1:k,joinline))
                scatter(X(setdiff(1:k,joinline),1),X(setdiff(1:k,joinline),2),'MarkerEdgeColor',icolor(i,:),'Marker',markerstr(i))
            end
        end
        
    case false % then input is the matrix of PNS coordinates
        
        [kk n]=size(data); k = (kk+1)/2 + 1; d = 2;
        PNS = varargin{2};
        
        if nargin < 3
            joinline = [1:k 1];
        else
            joinline = varargin{3};
        end
        icolor = repmat([0 0 0],n,1);
        markerstr = repmat('.',n,1);
        titlestring = 'shape configuration';
        
        if nargin > 3
            paramstruct = varargin{4};
            if isfield(paramstruct,'titlestring') ;    %  then change to input value
                titlestring = paramstruct.titlestring;
            end ;
            if isfield(paramstruct,'icolor') ;    %  then change to input value
                icolor = paramstruct.icolor;
                if length(icolor)==1 && icolor <  4
                    icolor = repmat('k',n,1);
                end
                if ischar(icolor) || size(icolor,1) == 1
                    icolor = repmat(icolor,n,1);
                end
                
            end ;
            if isfield(paramstruct,'markerstr') ;    %  then change to input value
                markerstr = paramstruct.markerstr;
                if ischar(markerstr)
                    markerstr = repmat(markerstr,n,1);
                end
            end ;
        end
        
        preshapeT= PNSe2s(data,PNS);
        H = Helmertsub(k);
        hold on;
        Xmat = zeros(k*n,2);
        for i=1:n;
            X = H'*reshape(preshapeT(:,i),(k-1),d);
            Xmat((1:k)+(i-1)*n,:) = X;
            plot(X(joinline,1),X(joinline,2),':','Color',icolor(i,:),'Marker',markerstr(i))
            if ~isempty(setdiff(1:k,joinline))
                scatter(X(setdiff(1:k,joinline),1),X(setdiff(1:k,joinline),2),'MarkerEdgeColor',icolor(i,:),'Marker',markerstr(i))
            end
        end
end
title(titlestring);

limtmp=axisSM(Xmat(:,1),Xmat(:,2));
axis equal
set(gca,'Xlim',limtmp(1:2),'Ylim',limtmp(3:4));

