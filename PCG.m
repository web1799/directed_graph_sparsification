function [x,iter,flag]=PCG(L_G,L_S,b,maxIter, tol)

x0=rand(length(b),1);
temp=L_G'*b;
r0=b-L_G*temp;
LsLs=L_S*L_S';
z0=lamgsolver(LsLs, r0);
p0=z0;
iter=1;
while iter<=maxIter
    temp1=L_G'*p0;
    alpha=r0'*z0/(temp1'*temp1);
    x=x0+alpha*p0;
    r1=r0-alpha*L_G*temp1;
    if norm(r1)<tol
        break;
    end
    z1=lamgsolver(LsLs, r1);
    beta=z1'*r1/(z0'*r0);   
    p0=z1+beta*p0;
    iter=iter+1;
    x0=x;
    r0=r1;
    z0=z1;
end
iter=iter-1;
if norm(r1)>tol
flag=1;
else
    flag=0;
end

end