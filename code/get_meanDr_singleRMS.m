function [meanDr,singleRMS]=get_meanDr_singleRMS(xyz)
global K K1 R F h P
xyz1=sqrt(sum(xyz.*xyz,2));
theta1=acos((xyz*K')./(K1*xyz1));
r1=(-2*P*cos(theta1)+sqrt(4*(P*cos(theta1)).^2+8*(sin(theta1)).^2*P*(R+h)))./(2*(sin(theta1)).^2);
point_d=sqrt(sum(xyz.*xyz,2));
Dr=r1-point_d;
meanDr=mean(Dr);%统计点到抛物面的平均距离
singleRMS=sqrt(mean(Dr.*Dr));