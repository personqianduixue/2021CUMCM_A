kkk=-30:5:30;
dh=kkk;
maxD_array=[0.3610,0.3834,0.4048,0.4272,0.4495,0.4731,0.4979,0.5228,0.5476,0.5725,0.5973,0.6272,0.6470];
figure('Position',[526.6,365,559.9,268.8]),plot(dh,maxD_array,'*-');hold on
bad_index=find(maxD_array>0.6);
scatter(dh(bad_index),maxD_array(bad_index),100,'rx');
legend('','不满足约束')
xlabel('h相对改变量/%');ylabel('主索点最大径向移动距离/m');beautiplot('small')
exportgraphics(gcf,'img\主索点最大径向移动距离随h改变量的变化.png','Resolution',400)