function basis = findBasisNormal(q)
% Return basis vectors for normal space at q. basis is a cell array of size
% n, one element per dimension of the curve

[n,T]=size(q);

qnorm=sqrt(sum(q.^2,1));
% guard against zero-velocity samples (q/|q| is taken as zero there)
qunit=q./qnorm;
qunit(:,qnorm==0)=0;

% remove the component along q using the same inner product as the
% tangent projection so the basis is exactly normal to q
qq=InnerProd_Q(q,q);
basis=cell(1,n);
for j=1:n
    h=repmat(q(j,:),n,1).*qunit;
    h(j,:)=h(j,:)+qnorm;
    basis{j}=h-q*InnerProd_Q(q,h)/qq;
end
