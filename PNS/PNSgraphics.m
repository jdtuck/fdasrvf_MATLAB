function [] = PNSgraphics(resmat,PNS,paramstruct)
% PNSgraphics :
% Inputs:
% resmat       : Matrix of residuals from PNSmain.m
%
% PNS          : a Matlab structure of parameters from PNSmain.m
%
% paramstruct - a Matlab structure of input parameters
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
%  fields       values
%
%  equalaxis    0 (default) - do NOT set axis equal length for scatterplot
%               1           - set axis equal length
%
%  titlestring  : string for title
%
%  !! The following is used in projplot1SM.m and projplot2SM.m.
%  !! This program will use these parameters only if there exist
%  !! the files, which are contained from the personal Statistical
%  !! Analyisis M-codes by J. S. Marron (University of North Carolina)
%  !! Otherwise, the program will use the built-in plot.m and hist.m.
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
%    idataconn        indices of data points to connect with line segments
%                     []  (default) for not connectiong any data points
%                     otherwise : x 2 matrix of pairs of indices of data points]
%                     (thus all intergers from 1,...,n).
%                     For time series data, this can give a clearer view of the
%                     data ordering by using [[1, 2];[2, 3];[3, 4];...].
%                     For bias adjustment, with matched pairs, each row should
%                     have the ind0ces of the matches.
%
%    idataconncolor   can be any of 'r', 'g', 'b', etc., to use that color for all
%                     default is 'k'
%                     or can be 2 for easy rainbow coloring,
%                         intended for time series of curves
%                         (Caution: this will use the first part of icolor,
%                          so might make most sense to use with icolor = 2,
%                          to avoid strange results)
%                     or can be color matrix, where the number of rows
%                     is the same as the number of rows of idataconn
%                         (has no effect for idataconn = [])
%
%    idataconntype    can be any of '-', '--', '-.', ':'
%                     default is '-'
%                     or can be character array (created using strvcat) of these,
%                     where the number of rows is the same as
%                     the number of rows of idataconn
%                         (has no effect for idataconn = [])
%
%    ibigdot          0  (default)  use Matlab default for dot sizes
%                     1  force large dot size in prints (useful since some
%                              postscript graphics leave dots too small)
%                              (Caution: shows up as small in Matlab view)
%                              Only has effect when markerstr = '.'
%
%    legendcellstr    cell array of strings for legend (nl of them),
%                     useful for (colored) classes, create this using
%                     cellstr, or {{string1 string2 ...}}
%                     Also a way to add a "title" to the first plot
%                             for this, use input of form:  {{string}}
%                     Also can indicate symbols, by just adding (at least
%                             for +,x.o) into the text
%
%    mlegendcolor     nl x 3 color matrix, corresponding to cell legends above
%                     (not needed when legendcellstr not specified)
%                     (defaults to black when not specified)
%
%    vaxlim        Vector of axis limits
%                        Use [] for default of all automatically chosen, by axisSM
%                        Use 1 for symmetrically chosen, by axisSM
%                            (often preferred for centered plots, as in PCA)
%                        Otherwise, must be 1 x 4 row vector of axis limits
%                        Note:  Use is generally not recommended,
%                        because defaults give "good visual impression
%                        of decomposition.  It is mostly intended for use
%                        with "very different" projections
%
%    iplotaxes        0 (default) do not plot axes
%                     1 plot axes as determined by direction vectors,
%                           using thin line type
%
%    See also

if nargin==0, help PNSgraphics, return, end     %   Display help

%  First set all parameters to defaults
%

base = 3 ;
equalaxis = 0;
titlestring = '';
%  Now update parameters as specified,
%  by parameter structure (if it is used)
%
if nargin > 2 ;   %  then paramstruct is an argument
    if isfield(paramstruct,'base') ;    %  then change to input value
        base = paramstruct.base;
    else
        base = 3; % spherical data
    end
    
    if isfield(paramstruct,'equalaxis') ;    %  then change to input value
        equalaxis = paramstruct.equalaxis;
    end ;
    
    if isfield(paramstruct,'titlestring') ;    %  then change to input value
        titlestring = paramstruct.titlestring;
    end ;
    
    if ~isfield(paramstruct,'icolor') ;    %  then set default value
        paramstruct.icolor = 2;
    end ;
    
    if ~isfield(paramstruct,'iplotaxes') ;    %  then set default value
        paramstruct.iplotaxes = 1;
    end ;
else
    paramstruct = struct(); % set for plots
end

% Initialize
[d n] = size(resmat);


% Proportion of variances explained
varPNS = sum(abs(resmat).^2, 2) / n;
cumvarPNS = cumsum(varPNS);
propcumPNS =cumvarPNS / cumvarPNS(end) * 100;
propPNS = varPNS / cumvarPNS(end) * 100;

% radii and relarive radii
radii = flipud(PNS.radii);
relradii = radii(1:end-1)./radii(2:end);

subplot(2,2,1)
plot(1:d, radii,':.');hold on;
plot(1:(d-1),relradii,':.r')
xlim([0.5,d+0.5]);xlabel('Dimension of PNS');ylabel('Radius of Nested Spheres')
ylim([0 1]);
legend('radius','rel. radius','Location','SouthEast')
title(titlestring);

subplot(2,2,2)
hbar = bar(1:d,propPNS,0.5);hold on;
plot(1:d, propcumPNS,'-r.');
set(hbar,'EdgeColor',[0.5 0.5 1])
xlim([0.5,d+0.5]);xlabel('Dimension of PNS');ylabel('Prop. of variance')

baseVec = {'ProcrustesMean' 'PNSmean' 'Geod. Mean' 'NO alignment'};
if base < 3 && isfield(PNS,'itype')
    switch PNS.itype
        case 'seq.test'
            title(['PNS, aligned to ' baseVec{base+1}]);
        case 'small'
            title(['PNS (no test), aligned to ' baseVec{base+1}]);
        case 'great'
            title(['PNG, aligned to ' baseVec{base+1}]);
    end
else
    title('PNS')
end

existSMfiles1 = exist('projplot1SM.m');
existSMfiles2 = exist('projplot2SM.m');

if d > 3
    subplot(2,3,4)
    if sqrt(varPNS(2)) > 1e-15
        if existSMfiles2 == 2
            projplot2SM(resmat(1:2,:),eye(2),paramstruct);
        else
            plot(resmat(1,:),resmat(2,:),'.')
        end
        xlabel('1st PNS');ylabel('2nd PNS');
        title(['var_1 = ' num2str(propPNS(1),'%3.2f') '%'])
        if equalaxis==1;xlimfix = get(gca,'XLim');axis equal;set(gca,'YLim',xlimfix);end
    else
        if existSMfiles1 == 2
            projplot1SM(resmat(1,:),1,paramstruct);
        else
            hist(resmat(1,:));ylabel('hist. of scores')
        end
        xlabel('1st PNS');
        title(['var_1 = ' num2str(propPNS(1),'%3.2f') '%'])
    end
    
    paramstruct.legendcellstr = {}; % no more legend to print
    paramstruct.mlegendcolor = [];
    
    subplot(2,3,5)
    if sqrt(varPNS(3)) > 1e-15
        if existSMfiles2 == 2
            projplot2SM(resmat([1,3],:),eye(2),paramstruct);
        else
            plot(resmat(1,:),resmat(3,:),'.')
        end
        xlabel('1st PNS');ylabel('3rd PNS');
        title(['var_2 = ' num2str(propPNS(2),'%3.2f') '%'])
        if equalaxis==1; axis equal;set(gca,'YLim',xlimfix,'XLim',xlimfix);end
    else
        if existSMfiles1 == 2
            projplot1SM(resmat(2,:),1,paramstruct);
        else
            hist(resmat(2,:));ylabel('hist. of scores')
        end
        xlabel('2nd PNS')
        title(['var_2 = ' num2str(propPNS(2),'%3.2f') '%'])
    end
    
    subplot(2,3,6)
    if sqrt(varPNS(3)) > 1e-15
        if existSMfiles2 == 2
            projplot2SM(resmat(2:3,:),eye(2),paramstruct);
        else
            plot(resmat(2,:),resmat(3,:),'.')
        end
        xlabel('2nd PNS');ylabel('3rd PNS');
        title(['var_3 = ' num2str(propPNS(3),'%3.2f') '%'])
        if equalaxis==1; axis equal;set(gca,'YLim',xlimfix,'XLim',xlimfix);end
    else
        title(['var_3 = ' num2str(propPNS(3),'%3.2f') '%'])
    end
else
    subplot(2,3,4:6)
    if sqrt(varPNS(2)) > 1e-15
        if existSMfiles2 == 2
            projplot2SM(resmat(1:2,:),eye(2),paramstruct);
        else
            plot(resmat(1,:),resmat(2,:),'.')
        end
        xlabel('1st PNS');ylabel('2nd PNS');
        title(['var_1 = ' num2str(propPNS(1),'%3.2f') '%'])
        if equalaxis==1;xlimfix = get(gca,'XLim');axis equal;set(gca,'YLim',xlimfix);end
    else
        if existSMfiles1 == 2
            projplot1SM(resmat(1,:),1,paramstruct);
        else
            hist(resmat(1,:));ylabel('hist. of scores')
        end
        xlabel('1st PNS')
        title(['var_1 = ' num2str(propPNS(1),'%3.2f') '%'])
    end
end