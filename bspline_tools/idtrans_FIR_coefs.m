function [b, c0] = idtrans_FIR_coefs(N)

% Use a precalculated table for low-order splines
idtrans_table = { {1,1}, {1,1}, {8, [1, 6, 1]}, {6, [1,4,1]}, ...
                {384, [1, 76, 230, 76, 1]}, {120, [1, 26, 66, 26, 1]},...
                {46080, [1, 722, 10543, 23548, 10543, 722, 1]},...
                {5040, [1, 120, 1191, 2416, 1191, 120, 1]} };

idtrans = idtrans_table{N +1};

b = idtrans{2};
c0 = idtrans{1};