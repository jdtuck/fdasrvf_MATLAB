function wproj = projectTangent(w,q,basis)

w=w-InnerProd_Q(w,q)*q;
% gram schmidt
b1=basis{1};
b2=basis{2};

basis1=b1/sqrt(InnerProd_Q(b1,b1));
b2=b2-InnerProd_Q(basis1,b2)*basis1;
basis2=b2/sqrt(InnerProd_Q(b2,b2));

bo{1}=basis1;
bo{2}=basis2;

wproj=w-InnerProd_Q(w,bo{1})*bo{1}-InnerProd_Q(w,bo{2})*bo{2};

