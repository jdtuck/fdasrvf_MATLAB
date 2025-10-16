function Y = Gram_Schmidt(X,varargin)
% GRAM_SCHMIDT orthonormalize basis
% -------------------------------------------------------------------------
% Orthonormalzie basis using Gram Schmidt
% 
% Usage: Y = Gram_Schmidt(X,varargin)
%
% Input:
% X: cell defining basis elements
% varargin{1} contains function name for Inner product
% varargin{2} contains optional parameter del
% varargin{3} contains optional parameter Phi
% varargin{4} contains optional parameter params.a
% varargin{5} contains optional parameter params.b
% 
% Output:
% Y: cell containing new basis

epsilon = 0.000005;
[~,cols] = size(X);
Y = size(X);
i = 1;
r = 1;
Y{1} = X{1};
while( i <= cols ) 
    tempvect = 0;
    for j = 1:i-1
        tempvect =  tempvect + feval(varargin{1},Y{j},X{r},varargin{3:end})*Y{j};    
    end
    Y{i} = X{r} - tempvect;
    tmp = feval(varargin{1},Y{i},Y{i},varargin{3:end});
    if( tmp > epsilon )
        Y{i} = Y{i}/sqrt(tmp);
        i = i + 1; 
        r = r + 1;
    else
        if(r < i)        
            r = r + 1;
        else
            break;
        end
    end
end
