%% 微信公众号：数学模型（ID: MATHmodels）
warning off
R = 300; D = 500;         % [m  ] 基准态时反射面为半径和口径
xp = 0; yp = 0; zp = -(1-0.466)*R;
xc = 0; yc = 0; zc = 0;

%% 读取附件1数据：主索节点的坐标(xm,ym,zm)和编号id
dat1 = importdata('附件1.csv');
id = dat1.textdata(2:end,1);
xyz = dat1.data;
xm = xyz(:,1); ym = xyz(:,2); zm = xyz(:,3);

%% 读取附件2数据：促动器下端点(地锚点)坐标、基准态时上端点(顶端)的坐标 xa, ya, za
dat2 = importdata('附件2.csv');
xyz = dat2.data;
xa = xyz(:,[1,4]); ya = xyz(:,[2,5]); za = xyz(:,[3,6]);

%% 读取附件3数据：反射面板对应的主索节点编号，即三角形顶点的 id
tris = table2array(readtable('附件3.csv'));
[~,tri] = ismember(tris,id);

%% 绘图

% 绘制反射面板
plot3(xm(tri(:,[1:3,1]))',ym(tri(:,[1:3,1]))',zm(tri(:,[1:3,1]))',...
    'color',[0.8,0.8,0.8]);
hold on

% 绘制下拉索和促动器
plot3([xm xa]', [ym ya]',[zm za]');

% 绘制口径圆
x = cosd(0:360); y = sind(0:360); z = ones(size(x));
H = -sqrt(R^2 -(D/2)^2);
plot3(x*D/2,y*D/2,z*H,'b')

% 绘制 C, P 两点
plot3(xc,yc,zc,'ro',xp,yp,zp,'mo')

axis equal
axis([300*[-1,1,-1,1],-350,50])
xlabel('x'); ylabel('y'); zlabel('z')

scatter3(select_point(:,1),select_point(:,2),select_point(:,3), 10,'filled')

