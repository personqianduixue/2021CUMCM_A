close all
A=cell2mat(Data1(:,2:4));
x=A(:,1);y=A(:,2);z=A(:,3);
scatter3(x,y,z)%散点图
axis equal
figure
[X,Y,Z]=griddata(x,y,z,linspace(min(x),max(x),1000)',linspace(min(y),max(y),1000),'v4');%插值
% pcolor(X,Y,Z);shading interp%伪彩色图
% axis equal
% figure,contourf(X,Y,Z) %等高线图
figure,surf(X,Y,Z);%三维曲面
axis equal
% figure,meshc(X,Y,Z)%剖面图
% axis equal
% view(0,0); 
figure,meshc(X,Y,Z);%s三维曲面（浅色）+等高线
axis equal
hidden off;