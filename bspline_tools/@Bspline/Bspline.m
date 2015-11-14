classdef Bspline
    
properties
    N           % spline order for each dimension
    c           % spline coefficient array
    bsh         % spline basis shift for each dimension
    btype       % spline basis type for each dimension, either 'basic' or 'dual'
    boundaryfun % boundary function for signal extension
end

methods
    function obj = Bspline(x, N, boundaryfun)
        if isa(x,'Bspline')
            obj.c = x.c;
            obj.N = x.N;
            obj.bsh = x.bsh;
            obj.btype = x.btype;
            obj.boundaryfun = x.boundaryfun;
        else
            ds = ndims(x);
     
            if(ds~=length(N) && length(N)>1)
                error('Bspline:ordererror', 'Spline order must be a scalar, or specified for each dimension of the input array');
            end
            
            if(nargin<3)
                boundaryfun = @mirrorbound_1;
            end

            % If a scalar order parameter is supplied, use it for all
            % dimensions (except singletons)
            if(length(N)<ds)
                Nsc = N;
                for dm = 1:ds;
                    if(size(x,dm)>1)
                        N(dm) = Nsc;
                    else
                        N(dm) = 0;
                    end
                end
            end

            for d = 1:ds
                % For 0,1 order the direct transform samples equals
                % the signal samples
                if(N(d)>1)
                    x = bsplineNdtrans(x, N(d), boundaryfun);
                end
                obj.btype{d} = 'basic';
                x = shiftdim(x, 1);
            end

            obj.c = x;
            obj.N = N;
            % Direct transform above gives symmetric basic spline
            obj.bsh = zeros(size(N));
            obj.boundaryfun = boundaryfun;

        end
    end
    
    function f = subsref(obj, sub)
        switch sub.type
            case '()'
                if not(length(sub.subs)==1 && isvector(obj.c)) && (length(sub.subs) ~= ndims(obj.c))
                    error('Bspline:subsref','Invalid number of subscripts')
                end
                
                % Subscripts are interpreted as input sample points. Thus real numbers are valid subscripts. 
                %
                % (y,x) = G(n,m)
                %
                
                 obj = obj.basis('basic');
                
                 f = obj.c;
                     
                 nds = length(sub.subs);
                 
                 % If all subscripts are vectors, use separable evaluation
                 % of splines
                 if(all(cellfun(@isvector,sub.subs)))
                     if(nds>1)
                         for dm = 1:nds
                            f = evalBsplineN1dim( f, obj.N(dm), sub.subs{dm}, obj.bsh(dm), obj.boundaryfun);
                            f = shiftdim(f, 1);
                         end
                     else
                         % evaluate 1st non-singleton dimension
                         dm = fnsdim(obj);
                         f = evalBsplineN1dim( f, obj.N(dm), sub.subs{1}, obj.bsh(dm), obj.boundaryfun);
                     end
                 else
                     f = evalBsplineNpoints( obj.c, obj.N, sub.subs, obj.bsh, obj.boundaryfun );                     
                 end
                
            case '.'
                if(strcmp(sub.subs,'c'))
                    f = obj.c;
                elseif(strcmp(sub.subs,'bsh'))
                    f = obj.bsh;
                elseif(strcmp(sub.subs,'N'))
                    f = obj.N;
                end
                
            otherwise
                error('Bspline:subsasgn',...
                      'Not a supported subscripted assignment')
        end
    end
    
    function sz = size(obj)
        sz = size(obj.c);
    end
    
    function y = double(obj)
        
        sub.type = '()';
        
        sz = size(obj.c);
        dms = 1:length(sz);
        
        for dm = dms
            sub.subs{dm} = 1:sz(dm);
        end
        
        y = subsref(obj, sub);
    end
    
    function dm = fnsdim(obj)
        % First non-signleton dimension
        dm = find(size(obj.c)>1, 1);
    end
    
    function nobj = diff(obj, ndif, dm)        
        if(nargin<2)
            ndif = 1;
        end
        
        if(nargin<3)
            dm = fnsdim(obj); % First non-singleton dimension
        end
        
        nobj = Bspline(obj);
        
        nN = obj.N;
        nN(dm) = nN(dm)-ndif;
        nobj.N =  nN;
        
        % New basis shift
        nbsh = obj.bsh;
        nbsh(dm) = mod(nbsh(dm)+0.5*ndif,1);
        nobj.bsh = nbsh;
        
        % Extend signal using boundary function
        nc = obj.boundaryfun(obj.c, (-ndif+1):(size(obj.c,dm)+ndif), dm);
        nc = diff(nc,ndif,dm);
        
        % Output subscripts
        dims = ndims(nc);
        out_subs = cell(1,dims);
        for dim = 1:dims
            out_subs{dim} = ':';
        end
        out_subs{dm} = 1:size(obj.c,dm);
        
        nobj.c = nc(out_subs{:});
    end
    
    function grd = gradient(obj)
        dms = ndims(obj.c);
        
        grd = cell([dms,1]);
        
        for dm = 1:dms
            grd{dm} = diff(obj, 1, dm);
        end
    end
    
    function igr = integral(obj, dm)

        if(nargin<2)
            dm = fnsdim(obj);
        end
        
        igr = Bspline(obj);

        % New spline order
        nN = obj.N;
        nN(dm) = nN(dm)+1;
        igr.N =  nN;        
        
        % New basis shift
        nbsh = obj.bsh;
        nbsh(dm) = rem(nbsh(dm)-0.5,1);
        igr.bsh = nbsh;        
        
        igr.c = cumsum(obj.c, dm);
    end
    
    % Convert between dual and basic b-spline basis
    function dobj = dual(obj, dm)
        
        dms = ndims(obj.c);
        
        if(nargin>1)
            % Find dual of a specific dimension
            dual_dim = zeros(1,dms);
            dual_dim(dm) = true;
        else
            % Find dual of all dimensions
            dual_dim = ones(1,dms);
        end       
        
        if(not(mod(obj.N,2)))
            error('Bspline:dual',...
                          'Duals are only supported for odd-order splines')
        end
        
        btype_new = cell(1,dms);
        
        for dm = 1:dms
            if(dual_dim(dm))
                if(strcmp('dual',obj.btype{dm}))
                    btype_new{dm} = 'basic';
                elseif(strcmp('basic',obj.btype{dm}))
                    btype_new{dm} = 'dual';
                else
                    error('Bspline:dual',...
                      'unknown basis type %s', obj.btype{dm})
                end
            else
                btype_new{dm} = obj.btype{dm};
            end
        end
        
        dobj = obj.basis(btype_new);
    end
    
    function eq = has_basis(obj, btype)
        dms1 = ndims(obj.c);
        dms2 = length(btype);
        
        if(dms1~=dms2)
            eq = false;
            return 
        end
        
        for dm = 1:dms1
            if(not(strcmp(obj.btype{dm},btype{dm})))
                eq = false;
                return
            end
        end
       
        eq = true;
    end
    
    function obj2 = basis(obj, btype, dm)
     
        dms = ndims(obj.c);
                    
        if(nargin>2)
            % Change basis of a specific dimension
            btype_new = obj.btype;
            btype_new{dm} = btype;
        else
            % Change basis of all dimensions
            if iscell(btype)
                btype_new = btype;
            else
                btype_new = cell(1,dms);
                for dm = 1:dms
                    btype_new{dm} = btype;
                end
            end
        end
        
        if(has_basis(obj,btype_new))
            % Basis unchanged
            obj2 = obj;
            return
        end
        
        obj2 = Bspline(obj);
        dc = obj2.c;
        
        for dm = 1:dms
            if(strcmp('dual',btype_new{dm}) && strcmp('basic',obj2.btype{dm}))
                dc = bsplineNidtrans(dc, obj2.N(dm)*2 + 1, obj.boundaryfun);
                obj2.btype{dm} = 'dual';
            elseif(strcmp('basic',btype_new{dm}) && strcmp('dual',obj2.btype{dm}) )
                dc = bsplineNdtrans(dc, obj2.N(dm)*2 + 1, obj.boundaryfun);
                obj2.btype{dm} = 'basic';
            end
            dc = shiftdim(dc, 1);
        end
        obj2.c = dc;
        
    end
    
    function robj = reduce(obj, dm)
        dms = ndims(obj.c);
        
        if(nargin>1)
            % Reduce a specific dimension
            reduce_dim = zeros(1,dms);
            reduce_dim(dm) = true;
        else
            % Reduce all non-singleton dimensions
            reduce_dim = size(obj.c)>1;
        end  
        
        robj = basis(obj, 'dual');
        cn = robj.c;
                    
        for dm = 1:dms
            if(reduce_dim(dm))
                % Prefilter
                [u2N, c0] = u2N_FIR_coefs(obj.N(dm));
                cn = 0.5*c0*filterFIR(u2N, cn, obj.boundaryfun);

                % Subsample
                subs = cell(1,dms);
                for d = dms
                    subs{d} = ':';
                end
                subs{1} = 1:2:size(cn,1);

                cn = cn(subs{:});
            end
            cn = shiftdim(cn,1);
        end
        
        robj.c = cn;
    end
    
    function eobj = expand(obj, dm)
        dms = ndims(obj.c);
        
        if(nargin>1)
            % Expand a specific dimension
            expand_dim = zeros(1,dms);
            expand_dim(dm) = true;
        else
            % Expand all non-singleton dimensions
            expand_dim = size(obj.c)>1;
        end
        
        eobj = basis(obj, 'basic');
        cn = eobj.c;
        
        for dm = 1:dms
            if(expand_dim(dm))
                % Find size and preallocate
                cm_size = size(cn);
                cm_size(1) = cm_size(1)*2;

                cm = zeros(cm_size);
                
                % Upsample
                subs = cell(1,dms);
                for d = dms
                    subs{d} = ':';
                end
                subs{1} = 1:2:cm_size(1);
                
                cm(subs{:}) = cn;
                
                % Postfilter
                [u2N, c0] = u2N_FIR_coefs(obj.N(dm));
                cn = c0*filterFIR(u2N, cm, obj.boundaryfun);
            end
            cn = shiftdim(cn, 1);
        end
        
        eobj.c = cn;
    end
    
end
    
end
