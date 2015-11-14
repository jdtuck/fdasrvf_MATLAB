classdef MscaleBspline
    
properties
    b_level
end

methods
    function obj = MscaleBspline(x_im, N_order, M_scaling, K_levels)
        obj.b_level = cell(1,K_levels);
        obj.b_level{1} = Bspline(x_im, N_order);
        
        for k = 2:K_levels
            obj.b_level{k} = reduce(obj.b_level{k-1}, M_scaling);
        end
    end
end

end