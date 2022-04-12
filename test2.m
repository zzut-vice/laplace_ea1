%% LAP测试2
close all
clear

input = stlread('D:\recent\model0914opt_use\grid0914V1.stl');
Points = input.Points;
ConnectivityList = input.ConnectivityList;
%% 检查连通性
%% 检查重复点
link_th = 1e-6;
temp2 = [];
for i = 1:length(input.Points)
    temp1 = sum(abs(input.Points - input.Points(i,:)),2);
    temp1 = temp1<link_th;
    temp1(i) = false;
    if any(temp1)
        disp([i find(temp1)'])
        temp2 = [temp2;i find(temp1)'];
    end
end

% 以较小索引替换较大的索引,暂仅限成对替换
for i = 1:size(temp2,2)
    idx_min = min(temp2(i,:));
    idx_max = max(temp2(i,:));
    ConnectivityList(ConnectivityList==idx_max) = idx_min;
    ConnectivityList(ConnectivityList>idx_max) = ConnectivityList(ConnectivityList>idx_max)-1;
    temp2(temp2>idx_max) = temp2(temp2>idx_max)-1;
    Points(idx_max,:) = [];
end

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
%% 构建del，X-axis
x_origin = Points(:,1);
del = L*x_origin;

%% 建立方程 L_x * x = del_x
idx_moving = randi(length(x_origin),50,1);% 锚点
x_moving = x_origin(idx_moving) + 100*rand(length(idx_moving),1);% 移动点位置
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

[y1]=sor(L_x'*L_x,L_x'*del_x,x0,1.8,1e1);% x = sor(A,b,x0,w,tol)
y2 = bicgstab(L_x'*L_x,L_x'*del_x,1e-6,10000);
% [y,n]=gauseidel(L_x'*L_x,L_x'*del_x,x0,1e-3); % (A,b,x0,ep)  
% y=GaussianSolver(L_x'*L_x,L_x'*del_x);
% [y,r] = linsolve(L_x'*L_x,L_x'*del_x);
%%
c= L_x'*L_x;
a = inv(c);


