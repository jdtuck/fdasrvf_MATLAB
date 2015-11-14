function test_suite = testBsplineN %#ok<STOUT>
  initTestSuite;
end


function testbspline0() %#ok<DEFNU>
    x = -4:1/12:4;
    
    y0 = bspline0(x);
    yn = bsplineGeneralN(x,0);
    
    
    assertVectorsAlmostEqual(y0,yn);
end

function testbspline1() %#ok<DEFNU>
    x = -4:1/12:4;
    y0 = bspline1(x);
    yn = bsplineGeneralN(x,1);
    
    assertVectorsAlmostEqual(y0,yn);
end

function testbspline2() %#ok<DEFNU>
    x = -4:1/12:4;
    y0 = bspline2(x);
    yn = bsplineGeneralN(x,2);
    
    assertVectorsAlmostEqual(y0,yn);
end

function testbspline3() %#ok<DEFNU>
    x = -4:0.05:4;
    
    y3 = bspline3(x);
    yn = bsplineGeneralN(x,3);
    
    
    assertVectorsAlmostEqual(y3,yn); 

end

function testbspline4() %#ok<DEFNU>
    x = -4:0.05:4;
    
    y4 = bspline4(x);
    yn = bsplineGeneralN(x,4);
    
    
    assertVectorsAlmostEqual(y4,yn); 

end

function testBspline3Kern() %#ok<DEFNU>
    k = -1:1;

    b = bsplineNkernel(k,3,1);
    bref = [1,4,1]/6;
    
    assertVectorsAlmostEqual(b,bref);
    
end

function testbspline3dtrans3vsNdtrans() %#ok<DEFNU>
    s = repmat([1,2,3,4,1,2,3,4],[1,20]);
    cN = bsplineNdtrans( s, 3 );
    c3 = bspline3dtrans( s );
    
    assertVectorsAlmostEqual(cN,c3);
end

function testBspline0dtrans() %#ok<DEFNU>
    % Trivial case, but it should not fail
     s = repmat([1,2,3,4,1,2,3,4],[1,20]);
     c = bsplineNdtrans(s,0);
     b = bsplineNkernel(0,0,1);
     srec = filterFIR(b,c);
     
     assertVectorsAlmostEqual(s,srec);
end

function testbspline3dtrans() %#ok<DEFNU>
     s = repmat([1,2,3,4,1,2,3,4],[1,20]);
     c = bspline3dtrans(s);
     b = bsplineNkernel(-1:1,3,1);
     srec = filterFIR(b,c);
     
     assertVectorsAlmostEqual(s,srec);
end

function testBspline4dtrans() %#ok<DEFNU>
     s = repmat([1,2,3,4,1,2,3,4],[1,20]);
     c = bsplineNdtrans( s, 4);
     b = bsplineNkernel(-2:2,4,1);
     srec = filterFIR(b,c);
     
     assertVectorsAlmostEqual(s,srec);
end

function testCalcBsplineRootTable() %#ok<DEFNU>
    %[ctest, ztest, btest] = calcBsplineRootTable( );
    zref = {[],[],-0.171573,-0.267949,[-0.361341,-0.0137254]',[-0.430575,-0.0430963]',...
            [-0.488295,-0.0816793, -0.00141415]', [-0.53528,-0.12255,-0.00914869]'};
   
    for tn = 1:length(zref)
        [btest, ctest] = idtrans_FIR_coefs(tn-1);
        
        ntail = (length(btest)-1)/2;
        bref = bsplineNkernel(-ntail:ntail,tn-1,1)*ctest;

        assertVectorsAlmostEqual(btest, bref);
        %assertVectorsAlmostEqual(ztest{tn}, zref{tn}, 'relative', (1e-5)*2, (1e-10)*4);
    end
end