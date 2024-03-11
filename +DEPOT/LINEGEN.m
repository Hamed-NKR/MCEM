n=5;

r=1:n;
%[cx,cy,cz]=meshgrid(dr,dr,dr);

cc=zeros(n,5);

for i=1:n
cc(i,:)=[i,1,r(i),0,0];
end

