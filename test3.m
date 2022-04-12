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





%% Lx=del
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
%% 构建del，3-axis
x_origin = Points;
del = L*x_origin;

%% 建立方程 L_x * x = del_x
idx_moving = randi(length(x_origin),50,1);% 锚点
moving_scale = 2000;
x_moving = x_origin(idx_moving,:) + (moving_scale*rand(length(idx_moving),3)-(moving_scale/2));% 移动点位置
% L_X
L_x = zeros(length(idx_moving),size(L,2)); 
for i = 1:length(idx_moving)
    L_x(i,idx_moving(i)) = 1;
end
L_x = [L;L_x]; 
% x0
x0 = [x_origin]; 

% del_x
del_x = [del;x_moving];

[y1]=sor(L_x'*L_x,L_x'*del_x,x0,1.8,1e2);% x = sor(A,b,x0,w,tol)
y2_temp = inv(L_x'*L_x);
y2 = y2_temp*L_x'*del_x;

% y2 = bicgstab(L_x'*L_x,L_x'*del_x,1e-6,10000);
% [y,n]=gauseidel(L_x'*L_x,L_x'*del_x,x0,1e-3); % (A,b,x0,ep)  
% y=GaussianSolver(L_x'*L_x,L_x'*del_x);
% [y,r] = linsolve(L_x'*L_x,L_x'*del_x);

trisurf(ConnectivityList,y1(:,1),y1(:,2),y1(:,3),'EdgeColor','none')
hold on
scatter3(x_moving(:,1),x_moving(:,2),x_moving(:,3),'*')
axis equal
%%
moving_scale = 1000;
x_moving = x_moving + (moving_scale*rand(length(idx_moving),3)-(moving_scale/2));% 移动点位置
del_x = [del;x_moving];
y2 = y2_temp*L_x'*del_x;
figure(1)
trisurf(ConnectivityList,y2(:,1),y2(:,2),y2(:,3),'EdgeColor','none')
hold on
scatter3(x_moving(:,1),x_moving(:,2),x_moving(:,3),'*')
axis equal
hold off


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
