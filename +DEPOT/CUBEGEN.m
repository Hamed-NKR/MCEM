n=3;

r=1:n;
%[cx,cy,cz]=meshgrid(dr,dr,dr);

cc=zeros(n^3,5);

for i=1:3
for j=1:3
for k=1:3
cc(i+(j-1)*n+(k-1)*n^2,:)=[i+(j-1)*n+(k-1)*n^2,1,r(i),r(j),r(k)];
end
end
end

