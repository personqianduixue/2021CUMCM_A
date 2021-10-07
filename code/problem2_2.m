clear,clc
load Data1.mat
load Data3.mat
global K K1 R F h P
Data3=string(Data3);
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
all_point_bianhao=string(Data1(:,1));
select_point_bianhao=all_point_bianhao(select_index);
theta=acos(costheta(select_index));
select_point_d=M1(select_index);
h=0.4380;
P=2*(F+h);
r=(-2*P*cos(theta)+sqrt(4*(P*cos(theta)).^2+8*(sin(theta)).^2*P*(R+h)))./(2*(sin(theta)).^2);
peak_point_xyz=(R+h)*K%抛物面顶点坐标
Delta_d0=r-select_point_d;%初始需要移动的距离Dalta_d0
m_arrow=select_point./(sqrt(sum(select_point.*select_point,2)));%主索点单位方向向量
select_point=select_point+m_arrow.*Delta_d0;
lianjietu=zeros(length(all_point_bianhao));%所有主索点连接图矩阵
Data3_weizhi=zeros(size(Data3));
for ii=1:4300
    for jj=1:3
        Data3_weizhi(ii,jj)=find(all_point_bianhao==Data3(ii,jj));
    end
end

select_weizhi=cell(1,length(select_index));
for ii=1:length(select_index)
    [row_index,col_index]=find(Data3_weizhi==select_index(ii));
    weizhi1=Data3_weizhi(row_index(col_index==1),:);
    weizhi2=Data3_weizhi(row_index(col_index==2),:);
    weizhi3=Data3_weizhi(row_index(col_index==3),:);
    weizhi2=[weizhi2(:,2),weizhi2(:,1),weizhi2(:,3)];
    weizhi3=[weizhi3(:,3),weizhi3(:,1),weizhi3(:,2)];
    weizhi=[weizhi1;weizhi2;weizhi3];
    willdelete=ones(1,size(weizhi,1));
    for iii=1:size(weizhi,1)
        for jjj=1:size(weizhi,2)
            if ~isempty(find(select_index==weizhi(iii,jjj)))
                weizhi(iii,jjj)=find(select_index==weizhi(iii,jjj));
            else
                willdelete(iii)=0;
            end
        end
    end
    weizhi=weizhi(find(willdelete==1),:);
    select_weizhi(ii)={weizhi};
end
%%

RMS_array=[10];
Sum_Delta_d=-Delta_d0;
for k=1:100
    All_RMS=[];
    Delta_d=[];
    for i=1:length(select_weizhi)
        now_weizhi=cell2mat(select_weizhi(i));
        all6_tonji_point_xyz=[];%每个主索点周围的统计点
        for ii=1:size(now_weizhi,1)
            each_weizhi=now_weizhi(ii,:);
            xyz1=select_point(each_weizhi(1),:);
            xyz2=select_point(each_weizhi(2),:);
            xyz3=select_point(each_weizhi(3),:);
            single_tonji_point_xyz=get_tonji_point(xyz1,xyz2,xyz3);%单个三角形内统计点坐标
            all6_tonji_point_xyz=[all6_tonji_point_xyz;single_tonji_point_xyz];
        end
        all6_2_tonji_point_xyz=[select_point(i,:);all6_tonji_point_xyz];%每个主索点周围的统计点+主索点本身
        [meanDr,singleRMS]=get_meanDr_singleRMS(all6_2_tonji_point_xyz);
        Delta_d=[Delta_d,meanDr];
        All_RMS=[All_RMS,singleRMS];
    end
    RMS_now=sqrt(mean(All_RMS.*All_RMS));
    RMS_array=[RMS_array,RMS_now];
    m_arrow=select_point./(sqrt(sum(select_point.*select_point,2)));%主索点单位方向向量
    select_point=select_point+m_arrow.*Delta_d';
    Sum_Delta_d=Sum_Delta_d-Delta_d';
    if RMS_array(end-1)-RMS_array(end)<0.001%收敛时为0.0146
        RMS_array(end)
        break
    end
end

%% 数据矫正
load('Data2.mat')
A_point=[0,0,0];
B_point=select_point;
D_point=source_point;
E_point=cell2mat(Data2(:,5:7));
E_point=E_point(select_index,:);
cos_theta=(sum(E_point.*D_point,2))./(sqrt(sum(E_point.*E_point,2)).*sqrt(sum(D_point.*D_point,2)));
b_length=sqrt(sum(B_point.*B_point,2));
DE=E_point-D_point;
a_length=sqrt(sum(DE.*DE,2));
c_length=b_length.*cos_theta+sqrt((b_length.*cos_theta).^2-b_length.^2+a_length.^2);
AE_length=sqrt(sum(E_point.*E_point,2));
Delta_D=AE_length-c_length;
chazhi=Delta_D-Sum_Delta_d;
max_D=max(abs(Sum_Delta_d))%0.4977
%% 作图
figure('Position',[307.4,401,560,245.6])
plot(1:3,RMS_array(2:4),'r*-')
xlabel('迭代次数'),ylabel('均方误差RSM变化');beautiplot
exportgraphics(gcf,'img\均方误差RSM随迭代次数的变化.png','Resolution',400)
%%
figure;scatter3(source_point(:,1),source_point(:,2),source_point(:,3),10);hold on
scatter3(select_point(:,1),select_point(:,2),select_point(:,3),10)
xlabel('X'),ylabel('Y');zlabel('Z')
scatter3(peak_point_xyz(1),peak_point_xyz(2),peak_point_xyz(3),40,'filled')%抛物面顶点
scatter3([0],[0],[0],40,'filled')
legend('原始节点','调整节点','抛物面顶点','球心')
axis equal
%% 保存数据
%save select_point.mat select_point
% resulttable=[{'X坐标（米）','Y坐标（米','Z坐标（米）'};num2cell(peak_point_xyz)];
% writecell(resulttable,'result.xlsx','Sheet','理想抛物面顶点坐标')
% result2=[select_point_bianhao,select_point];
% resulttable=[{'节点编号','X坐标（米）','Y坐标（米）','Z坐标（米）'};result2];
% writematrix(resulttable,'result.xlsx','Sheet','调整后主索节点编号及坐标')
% result3=[select_point_bianhao,Delta_D];
% resulttable=[{'对应主索节点编号','伸缩量（米）'};result3];
% writematrix(resulttable,'result.xlsx','Sheet','促动器顶端伸缩量')


%%
% for ii=1:4300
%     lianjietu(Data3_weizhi(ii,1),Data3_weizhi(ii,2))=1;
%     lianjietu(Data3_weizhi(ii,1),Data3_weizhi(ii,3))=1;
%     lianjietu(Data3_weizhi(ii,2),Data3_weizhi(ii,1))=1;
%     lianjietu(Data3_weizhi(ii,3),Data3_weizhi(ii,1))=1;
% end
% select_lianjietu=lianjietu(select_index,select_index);
% for i=1:length(all_point_bianhao)
%     
% end





