%灵敏性分析
clear,clc
h0=0.4380;
Last_RMS=[];
maxD_array=[];
bili_array=[];   
for kkk=30:30
    load Data1.mat
    load Data3.mat
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
    
    h=h0+kkk/100*h0;
    P=2*(F+h);
    r=(-2*P*cos(theta)+sqrt(4*(P*cos(theta)).^2+8*(sin(theta)).^2*P*(R+h)))./(2*(sin(theta)).^2);
    peak_point_xyz=(R+h)*K;%抛物面顶点坐标
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
        max_D=max(abs(Sum_Delta_d));
        if RMS_array(end-1)-RMS_array(end)<0.001%收敛时为0.0146
            Last_RMS=[Last_RMS,RMS_array(end)];
            max_D
            maxD_array=[maxD_array,max_D];
            break
        end
    end

%     old_point=select_point;
%     A=[cos(pi/2-alpha),-sin(pi/2-alpha),0;
%         sin(pi/2-alpha),cos(pi/2-alpha),0;
%         0,0,1];
%     B=[1,0,0;
%         0,cos(pi/2-beta),-sin(pi/2-beta);
%         0,sin(pi/2-beta),cos(pi/2-beta)];
%     C=B*A;%旋转矩阵
%     new_point=C*old_point';
%     new_point=new_point';
%     x=new_point(:,1);y=new_point(:,2);z=new_point(:,3);
%     [X,Y,Z]=griddata(x,y,z,linspace(min(x),max(x),1000)',linspace(min(y),max(y),1000),'v4');%插值
%     
%     size2cutslice = size(Z,1);
%     cutZ = zeros(1,size2cutslice);
%     cutX = zeros(1,size2cutslice);
%     cutY = zeros(1,size2cutslice);
%     row = floor(size2cutslice/2);
%     for i=1:size2cutslice
%         cutZ(i) = Z(row,i);
%         cutX(i) = X(row,i);
%         cutY(i) = Y(row,i);
%     end
%     
%     R=300.4;R1=-(1-0.466)*R;
%     DX=cutX(2:end)-cutX(1:end-1);
%     DZ=cutZ(2:end)-cutZ(1:end-1);
%     n_arraw=[-DZ;DX];%法向量
%     N_arraw=n_arraw./(sqrt(sum(n_arraw.*n_arraw,1)));%单位法向量
%     Rushe_arraw=[zeros(1,size(N_arraw,2));-ones(1,size(N_arraw,2))];%入射向量
%     Chushe_arraw=Rushe_arraw-2*sum(Rushe_arraw.*N_arraw).*N_arraw;%出射向量
%     K_=Chushe_arraw(1,:)./Chushe_arraw(2,:);%斜率倒数
%     jiao_X=K_.*(R1-cutZ(1:end-1))+cutX(1:end-1);%交点X坐标
%     bili1=sum(abs(jiao_X)<=0.5)/length(jiao_X);%在区域内比例
% 
%     size2cutslice = size(Z,1);
%     cutZ = zeros(1,size2cutslice);
%     cutX = zeros(1,size2cutslice);
%     cutY = zeros(1,size2cutslice);
%     row = floor(size2cutslice/2);
%     for i=1:size2cutslice
%         cutZ(i) = Z(i,row);
%         cutX(i) = X(i,row);
%         cutY(i) = Y(i,row);
%     end
% 
%     R=300.4;R1=-(1-0.466)*R;
%     DX=cutY(2:end)-cutY(1:end-1);
%     DZ=cutZ(2:end)-cutZ(1:end-1);
%     n_arraw=[-DZ;DX];%法向量
%     N_arraw=n_arraw./(sqrt(sum(n_arraw.*n_arraw,1)));%单位法向量
%     Rushe_arraw=[zeros(1,size(N_arraw,2));-ones(1,size(N_arraw,2))];%入射向量
%     Chushe_arraw=Rushe_arraw-2*sum(Rushe_arraw.*N_arraw).*N_arraw;%出射向量
%     K_=Chushe_arraw(1,:)./Chushe_arraw(2,:);%斜率倒数
%     jiao_X=K_.*(R1-cutZ(1:end-1))+cutY(1:end-1);%交点X坐标
%     bili2=sum(abs(jiao_X)<=0.5)/length(jiao_X);%在区域内比例 
%     bili=((bili1+bili2)/2)^2;
%     bili_array=[bili_array,bili];
end
%%
close all
kkk=-30:5:30;
dh=kkk;
figure('Position',[526.6,365,559.9,268.8]),plot(dh,maxD_array,'*-');hold on
bad_index=find(maxD_array>0.6);
scatter(dh(bad_index),maxD_array(bad_index),100,'rx');
legend('','不满足约束')
xlabel('h相对改变量/%');ylabel('主索点最大径向移动距离/m');beautiplot('small')
exportgraphics(gcf,'img\主索点最大径向移动距离随h改变量的变化.png','Resolution',400)
figure('Position',[526.6,365,559.9,268.8]),plot(dh,Last_RMS,'-*');
xlabel('h相对改变量/%');ylabel('迭代结束时的均方误差RSM/m');beautiplot('small')
exportgraphics(gcf,'img\迭代结束时的均方误差RSM随h改变量的变化.png','Resolution',400)
figure('Position',[526.6,365,559.9,268.8]),plot(dh,bili_array,'-*');
xlabel('h相对改变量/%');ylabel('接收比');beautiplot('small')
exportgraphics(gcf,'img\接收比随h改变量的变化.png','Resolution',400)
