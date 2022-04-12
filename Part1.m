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

%% Ls矩阵
A = sparse(length(Points),length(Points));% 连通矩阵
D = A;
for i = 1:length(ConnectivityList)
    A(ConnectivityList(i,:),ConnectivityList(i,:)) = 1;
end
D(logical(speye(length(A)))) = sum(A,2);
A(logical(speye(length(A)))) = 0;
L = speye(length(A)) - (inv(D)*A) ;
Ls = D-A;

%% 锚点
An = sparse(2,length(Points));
An(1,1) = 1;
An(2,10000) = 1;
%%
L_st = [Ls;An];
temp = L_st'*L_st;
I = speye(length(temp));

% temp1 = temp\I;
co = temp1*L_st';
% B1 = ;





