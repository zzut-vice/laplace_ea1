function [out] = CurCal(DeformedPoints)
%计算曲率
dim = [1 3];% 确定需计算的空间维度
% 选点
% 原方法
% curLenRange = [100];% 估计驻点区在X方向上的长度，应适当宽裕
% k = boundary(DeformedPoints(:,dim(1)),DeformedPoints(:,dim(2)),1);% 所选的二维面
% 
% temp =  DeformedPoints(k,1)>(max(DeformedPoints(:,1)) - curLenRange);
% k = k(temp);

curTopRange = 12;% 选取最靠前的X个点
k = boundary(DeformedPoints(:,dim(1)),DeformedPoints(:,dim(2)),1);% 所选的二维面
[~,curTopIdx] = sort(DeformedPoints(k,1),'descend');
temp =  curTopIdx(1:curTopRange);
k = k(temp);

% 转点云
ptCloud = zeros(length(k),3);
ptCloud(:,dim) = DeformedPoints(k,dim);
ptCloud = pointCloud(ptCloud);

ptCloudOut = pcdownsample(ptCloud, 'gridAverage', 1);%该行会影响顺序
while true
    
    k = boundary(ptCloudOut.Location(:,dim(1)),ptCloudOut.Location(:,dim(2)),1);%重新排序
    k = k(1:end-1);
    CPxy = ptCloudOut.Location(k,dim);

    % posi_arr = [ ];
    kappa_arr = [ ];
    % norm_arr = [ ];


    for num = 2:(length(CPxy)-1)
        x = CPxy(num-1:num+1,1);
        y = CPxy(num-1:num+1,2);
        [kappa,norm_l] = PJcurvature(x,y);
    %     posi_arr = [posi_arr;[x(2),y(2)]];
        kappa_arr = [kappa_arr;kappa];
    %     norm_arr = [norm_arr;norm_l];
    end
    kappa_arr = abs(kappa_arr);
    % 取平均
    cur_sigma = std(kappa_arr);% 计算标准差
    ave_idx = abs(kappa_arr-mean(kappa_arr))<(2*cur_sigma);% 剔除偏离标准差过多的曲率点
    out = mean(kappa_arr(ave_idx));
%     out = mean(kappa_arr); %原方法
    if isnan(out)
        ptCloudOut = ptCloud;
    else
        break
    end
end

% kappa_arr = smoothdata(kappa_arr,'gaussian',10); %不改变均值
%% 可视化
% quiver(posi_arr(:,1),posi_arr(:,2),...
%     kappa_arr.* norm_arr(:,1),kappa_arr.* norm_arr(:,2))
% hold on
% scatter(posi_arr(:,1),posi_arr(:,2))
% figure
% plot(kappa_arr)

end

function [kappa,norm_k] = PJcurvature(x,y)
    x = reshape(x,3,1);
    y = reshape(y,3,1);
    t_a = norm([x(2)-x(1),y(2)-y(1)]);
    t_b = norm([x(3)-x(2),y(3)-y(2)]);
    
    M =[[1, -t_a, t_a^2];
        [1, 0,    0    ];
        [1,  t_b, t_b^2]];

    a = M\x;
    b = M\y;

    kappa  = 2.*(a(3)*b(2)-b(3)*a(2)) / (a(2)^2.+b(2)^2.)^(1.5);
    norm_k =  [b(2),-a(2)]/sqrt(a(2)^2.+b(2)^2.);
end

