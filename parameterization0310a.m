% Laplacian Deformation 参数化
% 设计变量：简单抽几个点
% 设计约束：控制点不越过下述平面围成的空间
% 几何约束：1-尾端底面（X向）；2-对称面（Y向）；3-下底面（Z向）；
clear
%% 导入几何
input = stlread('grid_half_0914V1.stl');
Points = input.Points;
ConnectivityList = input.ConnectivityList;

%% 网格问题，需要去除部分重复点
link_th = 10; % 200 is ok~
temp2 = [];
idx_todelete = [];
Points_old = Points;
% 寻找允许容差内的重复点
for i = 1:length(Points_old)
    if ~ismember(i,idx_todelete)
        temp1 = sum(abs(Points_old - Points_old(i,:)),2);
        temp1 = temp1<link_th;
        temp1(i) = false;
        temp1(idx_todelete) = false;
        if any(temp1)
            temp2 = find(temp1)';
            disp([i temp2])% 替换索引
            idx_todelete = [idx_todelete temp2];

            % 对temp2执行替换
            for j = 1:length(temp2)
                ConnectivityList(ConnectivityList==temp2(j)) = i;
            end
        end
    end
end
disp(['替换点数量' num2str(length(idx_todelete))])

% 在P、C中删去无效的点序号
idx_todelete = unique(idx_todelete);
[Points,ConnectivityList] = LD_clear_mesh(Points,ConnectivityList,idx_todelete);
%% 根据拓扑信息构建L，（ L * x = del）
L = LD_formL(Points,ConnectivityList);
del = L*Points;
%% 根据锚点anchor（约束点和设计点）补充矩阵，需要记录两种锚点的索引
% 形成L_X，后续根据L_x'*L_x*x = L_x'*del_x来操作
% 选定锚点时即可生成L_X，直接获取方阵L_x_inv = (L_x'*L_x)^-1
% del_x项需要配合移动锚点的坐标，所以需后续生成

% a-固定锚点（几何约束）：几何约束：
%函数capture_flat：输入点云、基准面位置、面法向、容差
% 1-尾端底面（X向）
cons_idx1 = capture_flat(Points,[0 0 0],[1 0 0],1e-3);
% 2-对称面（Y向）
cons_idx2 = capture_flat(Points,[0 0 0],[0 1 0],1e-3);
% 3-下底面（Z向）
cons_idx3 = capture_flat(Points,[0 0 0],[0 0 1],1e-3);

% b-移动锚点（设计变量点）：
temp = [-51800	-440	2200
    -49688	-1136	3637
    -49295	-1585	982
    -43000	-2700	2890
    -42334	-1280	5209
    -34276	-2306	4524
    -33230	-2700	1800
    -27273	-858	5323
    -22450	-3259	1369
    -22400	-3130	3450
    -14871	-4720	-2139
    -13397	-1978	4839
    -8280	-894	5316
    -7561	-7017	2810
    -7315	-3429	3183
    -6884	-6747	1813
    -5361	-4086	1693
    -4008	-11351	2299
    -3217	-7097	2544
    -2099	-1510	5117];
moving_idx = [];
for i = 1:size(temp,1)
    moving_idx = [moving_idx;capture_flat(Points,temp(i,:))];
end


% 查看捕捉情况
if true
    trisurf(ConnectivityList,Points(:,1),Points(:,2),Points(:,3),'EdgeColor','none')
    hold on
    scatter3(Points(cons_idx1,1),Points(cons_idx1,2),Points(cons_idx1,3))
    scatter3(Points(cons_idx2,1),Points(cons_idx2,2),Points(cons_idx2,3))
    scatter3(Points(cons_idx3,1),Points(cons_idx3,2),Points(cons_idx3,3))
    scatter3(Points(moving_idx,1),Points(moving_idx,2),Points(moving_idx,3),'r*')
    legend(["几何体";"约束X";"约束Y";"约束Z";"设计点"])
end

% c-补充矩阵，构建L_x
% 好像需要对各个维度分别建立矩阵
% 1-X方向
temp_idx = [cons_idx1;moving_idx];
L_x1 = zeros(size(temp_idx,1),size(L,2)); 
for i = 1:length(temp_idx)
    L_x1(i,temp_idx(i)) = 1;
end
L_x1 = [L;L_x1];
L_x1_inv = full((L_x1'*L_x1)\speye(size(L_x1,2)));
del_x1_s = [del(:,1);Points(cons_idx1,1)];% 补充的是锚点的坐标，此处暂不完全，仅加入约束项
% 2-Y方向
temp_idx = [cons_idx2;moving_idx];
L_x2 = zeros(size(temp_idx,1),size(L,2)); 
for i = 1:length(temp_idx)
    L_x2(i,temp_idx(i)) = 1;
end
L_x2 = [L;L_x2];
L_x2_inv = full((L_x2'*L_x2)\speye(size(L_x2,2)));
del_x2_s = [del(:,2);Points(cons_idx2,2)];
% L_x2_inv = inv(L_x2'*L_x2);
% 3-Z方向
temp_idx = [cons_idx3;moving_idx];
L_x3 = zeros(size(temp_idx,1),size(L,2)); 
for i = 1:length(temp_idx)
    L_x3(i,temp_idx(i)) = 1;
end
L_x3 = [L;L_x3]; 
L_x3_inv = full((L_x3'*L_x3)\speye(size(L_x3,2)));
del_x3_s = [del(:,3);Points(cons_idx3,3)];

% 变形函数句柄
LD = @(x) LD_deformation(x,...
    L_x1_inv,L_x2_inv,L_x3_inv,...
    L_x1,L_x2,L_x3,...
    del_x1_s,del_x2_s,del_x3_s);
%% test，对移动锚点加入扰动，变形测试
moving_points0 = Points(moving_idx,:);
moving_scale = 1000;
test_moving = (moving_scale*rand(size(moving_points0,1),3)-(moving_scale/2));
moving_points = moving_points0;
%%
moving_points = moving_points + test_moving;% 移动点位置
Points_df = LD(moving_points);

figure(1)
trisurf(ConnectivityList,Points_df(:,1),Points_df(:,2),Points_df(:,3),'EdgeColor','none')
hold on
scatter3(moving_points(:,1),moving_points(:,2),moving_points(:,3),'r*')
axis equal
hold off
%% 优化函数
close all
addpath("aero\")
costfun = @(x) -aerofun(LD(moving_points0+reshape(x,[],3)),ConnectivityList);
moving_scale = 50;
options = optimoptions('fmincon','Display','iter');
x_opt = fmincon(costfun,zeros(1,length(moving_points0(:))),[],[],[],[],-moving_scale*ones(1,length(moving_points(:))),moving_scale*ones(1,length(moving_points(:))),[],options);
%% 显示计算结果
moving_points = moving_points0 + reshape(x_opt,[],3);
Points_df = LD(moving_points);
% 图1，变形示意
figure(1)
trisurf(ConnectivityList,Points_df(:,1),Points_df(:,2),Points_df(:,3),'EdgeColor','none')
hold on
scatter3(moving_points(:,1),moving_points(:,2),moving_points(:,3),'r*')
axis equal
hold off

% 图2，气动性能曲线
alpha_range = 0:25;
ma_range = [6 8 10];
[co_force_lift,co_force_drag] = CalPressure_alphama(Points,ConnectivityList,alpha_range,0,ma_range);
[co_force_lift_df,co_force_drag_df] = CalPressure_alphama(Points_df,ConnectivityList,alpha_range,0,ma_range);
temp1 = co_force_lift./co_force_drag;
temp2 = co_force_lift_df./co_force_drag_df;

figure(2)
for i = 1:length(ma_range)
subplot(length(ma_range),1,i)
plot(alpha_range,temp1(:,i),'r')
hold on
plot(alpha_range,temp2(:,i),'b')
grid on
xlabel('Alpha(°)')
ylabel('L/D')
legend(["Base";"Opt"])
title(['Ma = ' num2str(ma_range(i))])
end

figure(3)
mov_dis = sqrt(sum((Points_df - Points).^2,2));
pcshow(Points_df,mov_dis)
colormap(form_color('viridis',80))
c = colorbar;
title('Distance Between Two Point Clouds')

figure(4)
moving = reshape(x_opt,[],3);
trisurf(ConnectivityList,Points(:,1),Points(:,2),Points(:,3),'EdgeColor','none','FaceAlpha',0.4)
hold on
scatter3(Points(moving_idx,1),Points(moving_idx,2),Points(moving_idx,3),'r*')
quiver3(moving_points0(:,1),moving_points0(:,2),moving_points0(:,3),moving(:,1),moving(:,2),moving(:,3),0.3,'Color',[1 0 0])
axis equal

%%





function [Points,ConnectivityList] = LD_clear_mesh(Points,ConnectivityList,idx_todelete)
% 清理mesh，删去索引的P+重整C序列，无替换功能
Points_old = Points;
% 清理P
Points(idx_todelete,:) = [];
% 清理C
ConnectivityList((ConnectivityList(:,1)==ConnectivityList(:,2))|(ConnectivityList(:,1)==ConnectivityList(:,3))|(ConnectivityList(:,2)==ConnectivityList(:,3)),:) = [];
temp_skip = 0;
for i = 1:length(Points_old)
    if ~ismember(i,idx_todelete)
        ConnectivityList(ConnectivityList==i) = i-temp_skip;
    else
        temp_skip = temp_skip + 1;
    end
end
end

function L = LD_formL(Points,ConnectivityList)
% 根据拓扑信息构建L，（ L * x = del）
A = sparse(length(Points),length(Points));% 连通矩阵
D = A;
for i = 1:length(ConnectivityList)
    A(ConnectivityList(i,:),ConnectivityList(i,:)) = 1;
end
A(logical(speye(length(A)))) = 0;
D(logical(speye(length(A)))) = sum(A,2);
L = speye(length(A)) - (D\A) ;
disp(['构建L方阵规模：' num2str(length(L))])
end

function capture_idx = capture_flat(Points,flat_pos,flat_dir,th)
if nargin==2
    temp1 = sum(abs(Points - flat_pos),2);
    [~,capture_idx] = min(temp1);
elseif nargin==4
    if all(size(flat_dir) == [1 3])
    elseif all(size(flat_dir) == [3 1])
        flat_dir = flat_dir';
    else
        error('方向向量尺寸有问题')
    end
    flat_dir = flat_dir/norm(flat_dir);
    temp1 = Points - flat_pos;
    temp1 = abs(sum(temp1.*flat_dir,2));
    capture_idx = find(temp1<th);
end
end


