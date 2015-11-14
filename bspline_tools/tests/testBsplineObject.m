function test_suite = testBsplineObject %#ok<STOUT>
  initTestSuite;
end

function testBspline2Double() %#ok<DEFNU>
     s = repmat(magic(4),[4,4]);
     spobj = Bspline(s,3);
     
     stest = double(spobj);
     
     assertVectorsAlmostEqual(s, stest);
end

function testBsplineDual()
     s = repmat(magic(4),[10,10]);
     spobj = Bspline(s,3);
     spobj2 = dual(dual(spobj));
     assertVectorsAlmostEqual(spobj2.c, spobj.c);
end

function testBpslineDiff() %#ok<DEFNU>
    x = 1:101;
    Msub = 2;
    x_fine = 1:(1/Msub):101;
    
    
    % A trigonometric function with 1.order differential
    function y = trigfn(x)
        y = 100/(2*pi)*sin(1/100*2*pi*(x-1))-(x-1);
    end
    
    function y = trigfn_diff(x)
        y = cos(1/100*2*pi*(x-1))-1;
    end

    % Eval function on coarse grid and its derivative on fine grid
    s = trigfn(x);
    sd_fine = trigfn_diff(x_fine);
    
    tol = [0.1, 0.01, 1e-6, 1e-8, 1e-8, 1e-9];
    
    for N = 2:5
        
        % Fit spline to s with forward transform
        spl = Bspline(s,N);

        % Differentiate in spline domain
        spld = diff(spl);

        % Evaluate on fine grid
        spld_eval_fine = spld(x_fine);

        % Compare to analytically differentiated function on fine grid.
        % To avoid edge effects when checking accuracy, only central part
        % is used
        assertVectorsAlmostEqual(sd_fine(3*Msub*N:end-3*Msub*N), spld_eval_fine(3*Msub*N:end-3*Msub*N), 'relative', tol(N));
    end
    
end

function testBpslineIntegral() %#ok<DEFNU>
    x = 1:81;
    Msub = 2;
    x_fine = 1:(1/Msub):81;
    
    
    % A tap function with integral
    function y = fn_integral(x)
        y = 8*sqrt(pi)*erf(x/8 - 5)/2 ;
    end
    
    function y = fn(x)
        y = exp(-(x/8 - 5).^2);
    end

    % Eval function on coarse grid and its integral on fine grid
    s = fn(x);
    sint_fine = fn_integral(x_fine)-fn_integral(1);
    
    % NOTE: errors are mostly due to errors in approximation of fn,
    % not its integral
    tol = [1e-2, 1e-3, 1e-5, 1e-6, 1e-7, 1e-8];
    
    for N = 0:3
        
        % Fit spline to s with forward transform
        spl = Bspline(s,N);

        % integrate in spline domain
        splint = integral(spl);

        % Evaluate on fine grid
        splint_eval_fine = splint(x_fine)-splint(1);

        
        % NOTE: evaluating integral on higher-order splines will give
        % integral of mirrored parts as well. This should however just give
        % a constant offset for values far from the edge. For values close
        % to the edge things are a bit more muddy, should be looked into

        % Compare to analytically differentiated function on fine grid.
        % To avoid edge effects when checking accuracy, only central part
        % is used        

        assertVectorsAlmostEqual(sint_fine(3*Msub*(N+1):end-3*Msub*(N+1)), splint_eval_fine(3*Msub*(N+1):end-3*Msub*(N+1)), 'relative', tol(N+1));
    end
    
end

function testBsplineObj() %#ok<DEFNU>
    s = repmat([1,2,3,4,1,2,3,4],[1,20]);
    c = bsplineNdtrans( s, 4);
    spobj = Bspline(s,4);
    
    assertVectorsAlmostEqual(c, spobj.c);
end

function testBsplineObj3interp() %#ok<DEFNU>
    s = repmat([1,2,3,4,1,2,3,4],[1,20]);
    spobj = Bspline(s,3);
    
    stest = spobj(1:160);
    
    assertVectorsAlmostEqual(s, stest);
end

function testBsplineObj3interp2d() %#ok<DEFNU>
    s = repmat(magic(4),[5,4]);
    spobj = Bspline(s,3);
    
    stest = spobj(1:20,1:16);
    
    assertVectorsAlmostEqual(s, stest);
end

function testBsplineObj3interp2d2() %#ok<DEFNU>
    s = repmat(magic(4),[5,4]);
    spobj = Bspline(s,3);
    [x1,x2] = ndgrid(1:20,1:16);
    stest = spobj(x1,x2);
    
    assertVectorsAlmostEqual(s, stest);
end