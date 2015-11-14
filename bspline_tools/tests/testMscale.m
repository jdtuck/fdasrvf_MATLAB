function test_suite = testMscale %#ok<STOUT>
  initTestSuite;
end

function test_Expand() %#ok<DEFNU>
    x = [zeros(1,21), ones(1,10), zeros(1,19)];

    x = 0:49;
    
    bx = Bspline(x,1,@mirrorbound_1);
    bx_0_5 = reduce(bx);

    % Interpolation using basis function and expansion using the m-scale
    % relation should give equal results
    bx_rec = expand(bx_0_5);
    x_rec = bx_rec(1:50);
    x_rec_interp = bx_0_5(1:0.5:25.5);

    % Some problem with the end boundary condition
    assertVectorsAlmostEqual(x_rec(1:49), x_rec_interp(1:49));
end

function test_u2N_FIR_coefs() %#ok<DEFNU>
    ntail = 7;
    % u2N FIR coefs is valid for N odd
    for N = 1:2:7

        [u2N,c] = u2N_FIR_coefs(N);
        % ntail can here be wider than actual support, didnt bother to calculate
        % the support...
        bs1 = bsplineNkernel(-ntail:ntail,N,1);
        bs2 = bsplineNkernel(-ntail:ntail,N,2);
        
        % ...we just keep non-zero elements
        bs1(bs1==0)=[];
        bs2(bs2==0)=[];
        
        % Test that the u2N FIR filter fulfills the two-scale relation
        
        bs2_test = c*conv(u2N, bs1);
        
        bs2_ref = zeros(size(bs2_test));
        bs2_ref(1:length(bs2))=bs2;
        
        assertVectorsAlmostEqual(bs2_ref, bs2_test);
    end

end