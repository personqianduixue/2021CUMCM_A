%求解问题1理想抛物面h
clear,clc
load Data1.mat
R=300.4;
F=0.466*R;
alpha=0;beta=pi/2;
M=cell2mat(Data1(:,2:4));
K=[-cos(beta)*cos(alpha),-cos(beta)*sin(alpha),-sin(beta)];
M1=sqrt(sum(M.*M,2));
K1=sqrt(sum(K.*K));
costheta=(M*K')./(K1*M1);
select_index=find(costheta>cos(pi/6));
select_point=M(select_index,:);
theta=acos(costheta(select_index));
select_point_d=M1(select_index);
fitness=[];
for h=0:0.01:0.6
    P=2*(F+h);
    r=(-2*P*cos(theta)+sqrt(4*(P*cos(theta)).^2+8*(sin(theta)).^2*P*(R+h)))./(2*(sin(theta)).^2);
    r(1)=R+h;
    D=r-select_point_d;
    fitness=[fitness,sqrt(mean(D.*D))];
end
h=0:0.01:0.6;
figure('Position',[466.6,395.4,559.9,294.4]);
plot(h,fitness);xlabel('抛物面顶点与基准圆面距离h/m');ylabel('径向移动均方距离/m')
beautiplot('small')
exportgraphics(gcf,'img\问题1径向移动均方距离随h的变化.png','Resolution',400)
best_h=h(find(fitness==min(fitness)))%best_h =0.4400
%%
getfitness_p1(0.39,theta,select_point_d)%测试
getfitness_p1(0.44,theta,select_point_d)%测试
%%
f=@(h) getfitness_p1(h,theta,select_point_d);
[result,x]=huangjin(f,[0,0.6],100)%result = 0.1978;x =0.4444
%%
h=0.4444;
P=2*(F+h);
r=(-2*P*cos(theta)+sqrt(4*(P*cos(theta)).^2+8*(sin(theta)).^2*P*(R+h)))./(2*(sin(theta)).^2);
r(1)=R+h;
max(r.*sin(theta))
%%
function fit=getfitness_p1(h,theta,select_point_d) 
R=300.4;
F=0.466*R;
P=2*(F+h);
r=(-2*P*cos(theta)+sqrt(4*(P*cos(theta)).^2+8*(sin(theta)).^2*P*(R+h)))./(2*(sin(theta)).^2);
r(1)=R+h;
D=r-select_point_d;
fit=sqrt(mean(D.*D));
end

