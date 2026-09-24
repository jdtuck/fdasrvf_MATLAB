function wproj = projectTangent(w,q,basis)

w=w-InnerProd_Q(w,q)*q;
% gram schmidt
bo=cell(1,length(basis));
for i=1:length(basis)
    b=basis{i};
    for j=1:i-1
        b=b-InnerProd_Q(bo{j},b)*bo{j};
    end
    bo{i}=b/sqrt(InnerProd_Q(b,b));
end

wproj=w;
for i=1:length(bo)
    wproj=wproj-InnerProd_Q(w,bo{i})*bo{i};
end
