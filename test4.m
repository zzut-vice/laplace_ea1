%% LAP测试2
close all
clear

% input = stlread('D:\recent\model0914opt_use\grid0914V1.stl');
input = stlread('grid_half_0914V1.stl');

Points = input.Points;
ConnectivityList = input.ConnectivityList;

% 根据对称性删去部分

%% 检查连通性
%% 检查重复点
link_th = .1; % 200 is ok~
temp2 = [];
idx_todelete = [];
Points_old = Points;
for i = 1:length(Points_old)
    if ~ismember(i,idx_todelete)
        temp1 = sum(abs(Points_old - Points_old(i,:)),2);
        temp1 = temp1<link_th;
        temp1(i) = false;
        temp1(idx_todelete) = false;
        if any(temp1)
            temp2 = find(temp1)';
            disp([i temp2])
            idx_todelete = [idx_todelete temp2];

            % 对temp2执行替换
            for j = 1:length(temp2)
                ConnectivityList(ConnectivityList==temp2(j)) = i;
            end
        end
    end
end

% 在P、C中删去无效的点序号
idx_todelete = unique(idx_todelete);
[Points,ConnectivityList] = clear_mesh_idx(Points,ConnectivityList,idx_todelete);





%% L*x=del
% 构建L
A = sparse(length(Points),length(Points));% 连通矩阵
D = A;
for i = 1:length(ConnectivityList)
    A(ConnectivityList(i,:),ConnectivityList(i,:)) = 1;
end
A(logical(speye(length(A)))) = 0;
D(logical(speye(length(A)))) = sum(A,2);
L = speye(length(A)) - (D\A) ;
disp(length(L))
% 构建del，3-axis
x_origin = Points;
del = L*x_origin;

alpha = 10;% 约束强度，1-10
%% 锚点设置
mov_num = 50;% 移动锚点数量，自由度的根源
idx_mov = randperm(length(x_origin));
idx_mov = sort(idx_mov(1:mov_num),'ascend');% 定义移动锚点的索引
S = cal_S(idx_mov,size(x_origin,1));
anc_ini = x_origin(idx_mov,1);% 初始anchor 锚点

%% 加入随机扰动
moving_scale = 1000;
x_moving = x_origin(idx_mov,:) + (moving_scale*rand(length(idx_mov),3)-(moving_scale/2));% 移动点位置
% x_moving = x_moving + (moving_scale*rand(length(idx_mov),3)-(moving_scale/2));% 移动点位置

anc = x_moving(:,1);
%% 一旦确定锚点之后，可以生成loss(x)、loss_grad(x)句柄
loss_f = @(x) loss_function(L,x,del(:,1),alpha,S,anc);
loss_grad = @(x) loss_grad_function(L,x,del(:,1),alpha,S,anc);
loss_grad2 = @(x) loss_grad2_function(L,x,del(:,1),alpha,S,anc);

%
%%
fg2 = loss_grad2(x_origin(:,1));
fg2_inv = inv(fg2);
x_out = newton(loss_f,loss_grad,fg2_inv,x_origin(:,1),1e0,10000);
% x_out = gradient_descend(loss_f,loss_grad,x_origin(:,1),1e-2,1e0,10000);

function [Points,ConnectivityList] = clear_mesh_idx(Points,ConnectivityList,idx_todelete)
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

function S = cal_S(idx_moving,len)
S = sparse(len,len);
for i = 1:length(idx_moving)
    S(idx_moving(i),idx_moving(i)) = 1;
end
end


