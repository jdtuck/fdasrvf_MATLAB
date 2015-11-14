function test_suite = testPadding %#ok<STOUT>
  initTestSuite;
end

function testMirrorbound() %#ok<DEFNU>
    x = magic(4);
    
    ytest = mirrorbound_1(x,{-2:7,-2:7},1:2);
    
    yc = flipdim(flipdim(x,1),2);
    ylr = flipdim(x,2);
    ytb = flipdim(x,1);
    
    yref = [ yc(1:3,1:3), ytb(1:3,:), yc(1:3,2:4); ylr(:,1:3), x, ylr(:,2:4); yc(2:4,1:3), ytb(2:4,:), yc(2:4,2:4)];
    
    assertVectorsAlmostEqual(ytest,yref);
    
end