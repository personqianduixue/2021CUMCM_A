%求解问题2理想抛物面h
clear,clc
load Data1.mat
R=300.4;
F=0.466*R;
alpha=36.795/180*pi;beta=78.169/180*pi;
M=cell2mat(Data1(:,2:4));
K=[-cos(beta)*cos(alpha),-cos(beta)*sin(alpha),-sin(beta)];
M1=sqrt(sum(M.*M,2));
K1=sqrt(sum(K.*K));
costheta=(M*K')./(K1*M1);
select_index=find(costheta>cos(pi/6));
select_point=M(select_index,:);
source_point=select_point;
save source_point.mat source_point
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
plot(h,fitness)
best_h=h(find(fitness==min(fitness)))%best_h =0.4400
%%

f=@(h) getfitness_p1(h,theta,select_point_d);
[result,best_h2]=huangjin(f,[0,0.6],100)%result = 0.2032;best_h2 =0.4380
%%
% getfitness_p1(h,theta,select_point_d)
h=0.4380;
R=300.4;
F=0.466*R;
P=2*(F+h);
r=(-2*P*cos(theta)+sqrt(4*(P*cos(theta)).^2+8*(sin(theta)).^2*P*(R+h)))./(2*(sin(theta)).^2);
r(1)=R+h;
D=r-select_point_d;
%fit=sqrt(mean(D.*D));
mean(D)

function fit=getfitness_p1(h,theta,select_point_d) 
R=300.4;
F=0.466*R;
P=2*(F+h);
r=(-2*P*cos(theta)+sqrt(4*(P*cos(theta)).^2+8*(sin(theta)).^2*P*(R+h)))./(2*(sin(theta)).^2);
r(1)=R+h;
D=r-select_point_d;
fit=sqrt(mean(D.*D));
%fit=mean(D);
end

