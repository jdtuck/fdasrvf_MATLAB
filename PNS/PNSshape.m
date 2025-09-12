function [EuclidData, PNS] = PNSshape(data,paramstruct)
% PNSSHAPE : Principal Nested Spheres for planar SHAPE data.
%            - Draw graphics in current figure and returns the
%              Euclidean-type representation of the data with PNSmean,
%              and the PNS axes and radii
% example:
%   paramstruct = struct('base',0,'itype',0,...
%        'equalaxis',0,...
%        'icolor',2,'titlestring','title here');
%   [EuclidData, PNS] = PNSshape(data,paramstruct)
%
% Inputs:
% data    - k x m x n array of data, where k is the number of landmarks of
%               shape, m is the dimension of the space landmarks lie, and
%               n is the number of samples.
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
%  base      0 (default)'Procrustes' :  take the base w as 
%                                   the full Procrustes mean        
%
%            1  'PNSmean': take the base w as the PNSmean of the
%                          preshapes. 
%               
%            2  'gmean'  : take the base w (to which preshapes are rotated
%                          by an action in SO(2)) as geodesic mean of the
%                          preshapes.
%
%  itype     0  'seq.test' : (default) ordinary Principal Nested Sphere 
%                               with sequential tests.
%            1  'small'    : Principal Nested SMALL Sphere 
%            2  'great'    : Principal Nested GREAT Sphere (radius pi/2)
%  
%  alpha        0.05 (defualt) : size of Type I error allowed for each test
%               could be any number between 0 and 1.
%
%  R         100 (default) : number of bootsrap samples to be evaluated for
%            the sequential test.
%
%  equalaxis    0 (default) - do NOT set axis equal length for scatterplot
%               1           - set axis equal length
%
%  titlestring  : string for title
%  
%  printPvalue  0 (default) - do NOT print pvalues
%               1           - print pvalues
%
%  icolor           0  fully black and white version (everywhere)
%                     string (any of 'r', 'g', 'b', etc.) that single color
%                     1  (default)  color version (Matlab 7 color default)
%                     2  time series version (ordered spectrum of colors)
%                     nx3 color matrix:  a color label for each data point
%                             (to be used everywhere, except SiZer & QQ
%                              useful for comparing classes)
%
%  markerstr        Can be either a single string with symbol to use for marker,
%                         e.g. 'o' (default), '.', '+', 'x'
%                         (see "help plot" for a full list)
%                     Or a character array (n x 1), of these symbols,
%                         One for each data vector, created using:  strvcat
%
%  idataconn        indices of data points to connect with line segments
%                     []  (default) for not connectiong any data points
%                     otherwise : x 2 matrix of pairs of indices of data points]
%                     (thus all intergers from 1,...,n).
%                     For time series data, this can give a clearer view of the
%                     data ordering by using [[1, 2];[2, 3];[3, 4];...].
%                     For bias adjustment, with matched pairs, each row should
%                     have the ind0ces of the matches.
%
%  idataconncolor   can be any of 'r', 'g', 'b', etc., to use that color for all
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
%  idataconntype    can be any of '-', '--', '-.', ':'
%                     default is '-'
%                     or can be character array (created using strvcat) of these,
%                     where the number of rows is the same as
%                     the number of rows of idataconn
%                         (has no effect for idataconn = [])
%
%  ibigdot          0  (default)  use Matlab default for dot sizes
%                     1  force large dot size in prints (useful since some
%                              postscript graphics leave dots too small)
%                              (Caution: shows up as small in Matlab view)
%                              Only has effect when markerstr = '.'
%
%  legendcellstr    cell array of strings for legend (nl of them),
%                     useful for (colored) classes, create this using
%                     cellstr, or {{string1 string2 ...}}
%                     Also a way to add a "title" to the first plot
%                             for this, use input of form:  {{string}}
%                     Also can indicate symbols, by just adding (at least
%                             for +,x.o) into the text
%
%  mlegendcolor     nl x 3 color matrix, corresponding to cell legends above
%                     (not needed when legendcellstr not specified)
%                     (defaults to black when not specified)
%
%  vaxlim        Vector of axis limits
%                        Use [] for default of all automatically chosen, by axisSM
%                        Use 1 for symmetrically chosen, by axisSM
%                            (often preferred for centered plots, as in PCA)
%                        Otherwise, must be 1 x 4 row vector of axis limits
%                        Note:  Use is generally not recommended,
%                        because defaults give "good visual impression
%                        of decomposition.  It is mostly intended for use
%                        with "very different" projections
%
%  iplotaxes        0 (default) do not plot axes
%                     1 plot axes as determined by direction vectors,
%                           using thin line type
%
%  See also 

% Last updated Nov 6, 2009
% Sungkyu Jung

if nargin==0 && nargout==0, help PNSshape, return, end     %   Display help
%  First set all parameters to defaults
%
base = 0 ;
itype = 0 ;
alpha = 0.05;
R = 100;
%  Now update parameters as specified,
%  by parameter structure (if it is used)
if nargin > 1 ;   %  then paramstruct is an argument
    if isfield(paramstruct,'base') ;    %  then change to input value
        base = paramstruct.base;
    end ;
    
    if isfield(paramstruct,'itype') ;    %  then change to input value
        itype = paramstruct.itype ;
    end ;
    
    if isfield(paramstruct,'alpha') ;    %  then change to input value
        alpha = paramstruct.alpha;
    end ;
    if isfield(paramstruct,'R') ;    %  then change to input value
        R = paramstruct.R;
    end ;
else 
    paramstruct = struct(); % set for plots
end


baseVec = {'ProcrustesMean' 'PNSmean' 'Geod. Mean' 'axis'};
itypeVec = {'seq.test' 'small' 'great'};
disp(['Running PNSshape.m with base ' baseVec{base+1} ', itype = ' itypeVec{itype+1}]);

% k is the number of landmarks
% m is the dimension of the space where landmarks lie
% n is the size of sample
[k m n] = size(data);
% 
% if d~=2 ;
%     disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
%     disp('!!!   Error from PNSshape.m:      !!!') ;
%     disp('!!!   Shape is not 2-dimensional  !!!') ;
%     disp('!!!   Terminating execution       !!!') ;
%     disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
%     return ;
% end ;

% The preshape shape is 'dimpreshape'-dimensional sphere in ''dimpreshape'+1
% space
%dimpreshape = m*(k-1) - 1 ;

switch base
        case 0  % 'Procrustes'
        [rotatedData,ProcrustesMean] = procGPA(data);
        % Helmert submatrix for centering the shapes
        H = Helmertsub(k);
        % preshapes in matrix form
        preshapes = zeros(k-1,m,n);
        % preshapes in vector form
        Vpreshapes = zeros((k-1)*m,n);
        PmeanH = H*ProcrustesMean;
        PmeanZ = PmeanH/norm(PmeanH,'fro'); % Procrustes mean preshape
        for i=1:n;
            X = rotatedData(:,:,i);
            XH = H*X; % translation
            Z = XH/norm(XH,'fro'); % scaling
            
            % Now align
            rotZ=optrotPreshape(Z,PmeanZ);
            preshapes(:,:,i) = rotZ;
            Vpreshapes(:,i) = rotZ(:);
        end
        [EuclidData PNS]= PNSmain(Vpreshapes,itype,alpha,R);
        PNS.alignbase = PmeanZ;

    case 1 % 'PNSmean' % base w as the PNSmean
        [rotatedData,ProcrustesMean] = procGPA(data);
        % Helmert submatrix for centering the shapes
        H = Helmertsub(k);
        % preshapes in matrix form
        preshapes = zeros(k-1,m,n);
        % preshapes in vector form
        Vpreshapes = zeros((k-1)*m,n);
        PmeanH = H*ProcrustesMean;
        PmeanZ = PmeanH/norm(PmeanH,'fro'); % Procrustes mean preshape
        for i=1:n;
            X = rotatedData(:,:,i);
            XH = H*X; % translation
            Z = XH/norm(XH,'fro'); % scaling
            
            % Now align
            rotZ=optrotPreshape(Z,PmeanZ);
            preshapes(:,:,i) = rotZ;
            Vpreshapes(:,i) = rotZ(:);
        end
        [EuclidData PNS]= PNSmain(Vpreshapes,itype,alpha,R);
        PNSmean = PNS.mean;
        
        % Now rotate to the PNS mean until convergence is made

        diff = 1;cnt = 1;
        while diff > 1e-6
            for i=1:n
                zi=optrotPreshape(preshapes(:,:,i),reshape(PNSmean,k-1,m));
                preshapes(:,:,i) = zi;
                Vpreshapes(:,i) = zi(:);
            end
            [EuclidData PNS]= PNSmain(Vpreshapes,itype,alpha,R);
            diff = norm(PNS.mean - PNSmean,'fro');
            PNSmean = PNS.mean;
            disp(num2str(cnt));cnt = cnt+1;
            if cnt >= 10;
                disp(['iteration reached 10th step with diff ' num2str(diff)]);
                break;
            end
        end
        PNS.alignbase = PNSmean;
        

    case 2 % 'gmean'  base w as the geodesic mean of preshapes
        % Helmert submatrix for centering the shapes
        H = Helmertsub(k);
        % preshapes in matrix form
        preshapes = zeros(k-1,m,n);
        % preshapes in vector form
        Vpreshapes = zeros((k-1)*m,n);
        for i=1:n;
            X = data(:,:,i);
            XH = H*X; % translation
            Z = XH/norm(XH,'fro'); % scaling
            preshapes(:,:,i) = Z;
            Vpreshapes(:,i) = Z(:);
        end
        
        diff =1;cnt = 1;
        % initial geodesic mean
        mpreshape = reshape(geodmeanSk(Vpreshapes),k-1,m);
        
        while diff > 1e-10
            for i=1:n
                % rotate preshapes near the geodesic mean
                zi=optrotPreshape(preshapes(:,:,i),mpreshape);
                % undate rotated preshapes
                preshapes(:,:,i) = zi;
                Vpreshapes(:,i) = zi(:);
            end
            % updata geodesic mean based on rotated preshapes
            mpreshapeRot = reshape(geodmeanSk(Vpreshapes),k-1,m);
            % iterate until no further change in geodesic mean
            diff = norm(mpreshape - mpreshapeRot,'fro') ;           
            mpreshape = mpreshapeRot;
            cnt = cnt+1;
            if cnt > 10;
                break;
            end
        end
        [EuclidData PNS]= PNSmain(Vpreshapes,itype,alpha,R);
        PNS.alignbase = mpreshape;

   otherwise
        disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
        disp('!!! Error from PNSshape.m :                      !!!!')
        disp('!!! "base" does not match any option             !!!!')
        disp('!!!!Try 0,1, or 2                                !!!!')
        disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
        return;
end

% then draw plots

PNSgraphics(EuclidData,PNS,paramstruct)
