
 
% generate the data in domain [-1,1] 
N_points = 500;
X = linspace(-1,1,500);
Y = linspace(-1,1,500);
[meshX,meshY] = meshgrid(X,Y);
 
% construct the mesh solution by meshgrid of XY
meshU = sin(meshX.^2+meshY.^2);
 
% plot the mesh solution
figure('name','meshU')
mesh(meshX, meshY, meshU);
hold on
xlabel('x')
ylabel('y')
zlabel('u')
 
% extracting the diagonal parts from [-1,-1] to [1,1]
size2cutslice = size(meshU,1);
utrue_cut = zeros(size2cutslice, 1);
cutX_coord = zeros(size2cutslice, 1);
cutY_coord = zeros(size2cutslice, 1);
row = 30;
for i=1:size2cutslice
    meshu_cut(i,1) = meshU(i,i);
    cutX_coord(i,1) = meshX(i,i);
    cutY_coord(i,1) = meshY(i,i);
end
plot3(cutX_coord, cutY_coord, meshu_cut, 'r*', 'linewidth', 2)
hold on

figure;plot(cutX_coord,meshu_cut)
 
% extracting the cut line for a given xcoord
size2cutslice2 = size(meshU,1);
utrue_cut1 = zeros(size2cutslice2, 1);
cutX_coord1 = zeros(size2cutslice2, 1);
cutY_coord1 = zeros(size2cutslice2, 1);
row = 250;
for i=1:size2cutslice
    meshu_cut1(i,1) = meshU(i,row);
    cutX_coord1(i,1) = meshX(i,row);
    cutY_coord1(i,1) = meshY(i,row);
end
 
plot3(cutX_coord1, cutY_coord1, meshu_cut1, 'b+', 'linewidth', 2)
hold on
figure('name','cutline')
plot(meshu_cut, 'b--', 'linewidth', 2)
hold on
plot(meshu_cut1, 'm:', 'linewidth', 2)
hold on
 

% 对于沿固定的y轴方向，类似上述操作。
% 对于截取平行于xoy平面的曲线使用 matlab自带命令 contour 就好