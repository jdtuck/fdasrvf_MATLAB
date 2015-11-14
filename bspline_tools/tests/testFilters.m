function test_suite = testFilters %#ok<STOUT>
  initTestSuite;
end


function testFilterFIRmirror() %#ok<DEFNU>
     s = repmat([1,2,3,4,1,2,3,4],[1,5]);
     c = bspline3dtrans(s);
     b = bsplineNkernel(-1:1,3,1);
     srec1 = filterFIR(b,c);
     srec2 = conv([c(2) c c(end-1)],b,'valid');
     
     assertVectorsAlmostEqual(srec1,srec2);
end


function testSymm2ordAp() %#ok<DEFNU>
    s = repmat([1,2,3,4,1,2,3,4],[1,10]);
    
    z1 = -2 + sqrt(3);
    c1 = 6;
    
    y1 = c1*symm2ordAp_simple(s, z1);
    y2 = c1*filterAPsym2ord(s, z1);
    
    assertVectorsAlmostEqual(y1, y2);
end