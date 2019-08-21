function val = innerProdL2(f1,f2,t)

val=trapz(t,sum(f1.*f2,1));