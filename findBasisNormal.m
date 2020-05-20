function basis = findBasisNormal(q)
% Return basis vectors for normal space at q. basis is a cell array of size

[n,T]=size(q);

f1=zeros(n,T);
f2=zeros(n,T);
for i=1:T
    f1(:,i)=q(1,i)*q(:,i)/norm(q(:,i))+[norm(q(:,i));0];
    f2(:,i)=q(2,i)*q(:,i)/norm(q(:,i))+[0;norm(q(:,i))];
end
h3=f1;
h4=f2;
integrandb3=zeros(1,T);
integrandb4=zeros(1,T);
for i=1:T
    integrandb3(i)=q(:,i)'*h3(:,i);
    integrandb4(i)=q(:,i)'*h4(:,i);
end
b3=h3-q*trapz(linspace(0,1,T),integrandb3);
b4=h4-q*trapz(linspace(0,1,T),integrandb4);

basis{1}=b3;
basis{2}=b4;