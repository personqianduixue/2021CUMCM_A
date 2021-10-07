% load select_point.mat
% old_point=select_point;
load source_point.mat
old_point=source_point;
alpha=36.795/180*pi;beta=78.169/180*pi;
A=[cos(pi/2-alpha),-sin(pi/2-alpha),0;
    sin(pi/2-alpha),cos(pi/2-alpha),0;
    0,0,1];
B=[1,0,0;
    0,cos(pi/2-beta),-sin(pi/2-beta);
    0,sin(pi/2-beta),cos(pi/2-beta)];
C=B*A;%旋转矩阵
new_point=C*old_point';
new_point=new_point';
x=new_point(:,1);y=new_point(:,2);z=new_point(:,3);
figure;scatter3(old_point(:,1),old_point(:,2),old_point(:,3))
figure;scatter3(x,y,z)
%%
[X,Y,Z]=griddata(x,y,z,linspace(min(x),max(x),1000)',linspace(min(y),max(y),1000),'v4');%插值
%figure;meshc(X,Y,Z)
size2cutslice = size(Z,1);
cutZ = zeros(1,size2cutslice);
cutX = zeros(1,size2cutslice);
cutY = zeros(1,size2cutslice);
row = floor(size2cutslice/2);
for i=1:size2cutslice
    cutZ(i) = Z(row,i);
    cutX(i) = X(row,i);
    cutY(i) = Y(row,i);
end
figure;scatter(cutX,cutZ,5),xlabel('X'),ylabel('Z');title('ZX平面');magnify
% xlim(-150,150),ylim(-250,-310);axis equal

R=300.4;R1=-(1-0.466)*R;
DX=cutX(2:end)-cutX(1:end-1);
DZ=cutZ(2:end)-cutZ(1:end-1);
n_arraw=[-DZ;DX];%法向量
N_arraw=n_arraw./(sqrt(sum(n_arraw.*n_arraw,1)));%单位法向量
Rushe_arraw=[zeros(1,size(N_arraw,2));-ones(1,size(N_arraw,2))];%入射向量
Chushe_arraw=Rushe_arraw-2*sum(Rushe_arraw.*N_arraw).*N_arraw;%出射向量
K_=Chushe_arraw(1,:)./Chushe_arraw(2,:);%斜率倒数
jiao_X=K_.*(R1-cutZ(1:end-1))+cutX(1:end-1);%交点X坐标
bili1=sum(abs(jiao_X)<=0.5)/length(jiao_X)%在区域内比例 :0.9209

size2cutslice = size(Z,1);
cutZ = zeros(1,size2cutslice);
cutX = zeros(1,size2cutslice);
cutY = zeros(1,size2cutslice);
row = floor(size2cutslice/2);
for i=1:size2cutslice
    cutZ(i) = Z(i,row);
    cutX(i) = X(i,row);
    cutY(i) = Y(i,row);
end
figure;scatter(cutY,cutZ,5),xlabel('Y'),ylabel('Z');title('ZY平面');
% xlim(-150,150),ylim(-250,-310);axis equal
%%
R=300.4;R1=-(1-0.466)*R;
DX=cutY(2:end)-cutY(1:end-1);
DZ=cutZ(2:end)-cutZ(1:end-1);
n_arraw=[-DZ;DX];%法向量
N_arraw=n_arraw./(sqrt(sum(n_arraw.*n_arraw,1)));%单位法向量
Rushe_arraw=[zeros(1,size(N_arraw,2));-ones(1,size(N_arraw,2))];%入射向量
Chushe_arraw=Rushe_arraw-2*sum(Rushe_arraw.*N_arraw).*N_arraw;%出射向量
K_=Chushe_arraw(1,:)./Chushe_arraw(2,:);%斜率倒数
jiao_X=K_.*(R1-cutZ(1:end-1))+cutY(1:end-1);%交点X坐标
bili2=sum(abs(jiao_X)<=0.5)/length(jiao_X)%在区域内比例 0.9099
bili=((bili1+bili2)/2)^2%0.8380