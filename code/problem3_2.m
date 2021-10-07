R=300;
x=[];
for x0=1:0.001:150
    y0=-sqrt(R^2-x0^2);
    n=[-x0,-y0];
    e=[0,-1];
    n1=n/norm(n);
    en=sum(e.*n1);
    c=e-2*en*n1;
    x=[x,c(1)/c(2)*(-0.534*R-y0)+x0];
end
sum(abs(x)<0.5)/length(x)    