clc,clear
x=[];
R=300.4;
for x0=-150:0.1:150
    y0=-sqrt(R^2-x0.^2);
    n=[-x0;-y0];
    e=[0;-1];
    n1=n/norm(n);
    en=sum(e.*n1);
    c=e-2.*en.*n1;
    x=[x,c(1)/c(2).*(-0.534*R-y0)+x0];
end
x0=-150:0.1:150;
plot(x0,x)   
hold on
y11=[0.5;-0.5].*ones(2,length(x0));
plot(x0,y11,'r--')
xlabel('入射点横坐标');ylabel('交点横坐标')
beautiplot