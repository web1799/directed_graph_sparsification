function [x, GaussItr]=Gauss_Seidel(A, b, x,tol,maxiter)
plotGauss=[];
[n,m]=size(A);
GaussItr=0;
normVal=Inf; 
xnorm=norm(x);
while ((normVal>tol)&(GaussItr < maxiter))
    x_old=x;
%     x=x/norm(x)*xnorm;
%     x_old=x_old/norm(x_old)*xnorm;
    for i=1:n
        
        sigma=0;
        
        for j=1:i-1
                sigma=sigma+A(i,j)*x(j);
        end
        
        for j=i+1:n
                sigma=sigma+A(i,j)*x_old(j);
        end  
        x(i)=(1/A(i,i))*(b(i)-sigma);
    end
    
    GaussItr=GaussItr+1;
    normVal=norm(x_old-x);
    plotGauss=[plotGauss;normVal];
end
x=x/norm(x);
% x=x/norm(x)*xnorm;
% fprintf('Solution of the system is : \n%f\n%f\n%f\n%f\n%f in %d iterations',x,GaussItr);
end
