function [co_force_lift,co_force_drag,...
    co_force_axis,co_force_lateral,co_force_normal,...
    co_force_moment_roll,co_force_moment_pitch,co_force_moment_yaw,aero_pressure_center_ba,...
    S_r]...
    = CalPressure_alphama(Points,Connectivity,AOARange,beta,MaRange)

%% 初始化
ga = 1.4;
%% 外形
%% 构建面元
p0 = zeros(size(Connectivity));
n1 = zeros(size(Connectivity));
s = zeros(length(n1),1);
for i = 1:size(Connectivity,1)
    p0(i,:) = (Points(Connectivity(i,1),:) + Points(Connectivity(i,2),:) ...
        + Points(Connectivity(i,3),:)) /3;%各个网格的质心
    temp1 = cross(Points(Connectivity(i,3),:) - Points(Connectivity(i,1),:), ...
        Points(Connectivity(i,2),:) - Points(Connectivity(i,1),:));%叉乘方法确定外法向
    n1(i,:) = temp1/norm(temp1,2); %单位化的外法向
end
s = tri_area_mat(Connectivity,Points); % 计算三角网格面积（向量化）

S_r = sum(s)/1e6; % 参考面积 [m2]
L_r = (max(Points(:,1)) - min(Points(:,1)))/1e3;% 参考长度[m]
mass_center = sum(p0.*s)/sum(s); %形心替代质心


%% 容积率模块，仅适用于凸包
% VolumeEff = [];
% HeightGrid = abs(sum(n1.* (p0-InnerPoints),2));%三棱锥的高
% HeightGrid(PointsIndeformed.List) = 0;
% VolumeGrid = s .* HeightGrid./3;%三棱锥体积：a*h/3
% VolumeEff=6 * sqrt(pi) * sum(VolumeGrid) /(sum(s)^1.5);%容积率
%% 计算头部曲率，曲率越小越钝

%% 计算不同攻角下的受力
aero_pressure_center_ba = zeros(length(AOARange),length(MaRange));
co_force_axis = zeros(length(AOARange),length(MaRange)); % 轴向力系数
co_force_lateral = zeros(length(AOARange),length(MaRange)); % 侧向力系数
co_force_normal = zeros(length(AOARange),length(MaRange)); % 法向力系数
co_force_lift = zeros(length(AOARange),length(MaRange));% 升力系数
co_force_drag = zeros(length(AOARange),length(MaRange));% 阻力系数
co_force_moment_roll = zeros(length(AOARange),length(MaRange)); % 滚转系数
co_force_moment_pitch = zeros(length(AOARange),length(MaRange));% 俯仰系数
co_force_moment_yaw = zeros(length(AOARange),length(MaRange));% 偏航系数

for i = 1:length(AOARange)
    for j = 1:length(MaRange)
        Ma = MaRange(j);
        % 气流方向，无角度参数时为[-1 0 0]
        alpha = AOARange(i);
        n_stream = [-1 0 0]';
        n_stream = roty(alpha)*n_stream;
        n_stream = rotz(-beta)*n_stream;
        % 计算面元与气流夹角
        sub_angle_cos = n1*-n_stream;
        sub_angle = acos(n1*-n_stream);
        
        d = (2/(ga*(Ma^2)))*((((ga+1)^2*Ma^2)/(4*ga*Ma^2-2*(ga-1)))^(ga/(ga-1))*...
            ((1-ga+2*ga*Ma^2)/(ga+1))-1); % 驻点压力系数公式
        cp = zeros(length(s),1);
        % 压力系数计算：
        % 循环版本
%         for e_id = 1:length(s)
%             if s(e_id,1) == 0
%                 cp(e_id) = 0;
%             elseif sub_angle_cos(e_id) == 0
%                 % 与来流平行
%                 cp(e_id) = 0;
%             elseif sub_angle_cos(e_id) >0
%                 %迎风面
%                 cp(e_id) = d * sub_angle_cos(e_id) ^2;
%             else
%                 %背风面: Prandtl-Meyer method
%                 del = -(sub_angle(e_id) - pi/2);
%                 cp(e_id) = -(del.^2)*(ga/2 + 1/2).*((16./(del.^2*Ma^2*(ga + 1)^2) + 1)^(1/2) - 1);
%             end
%         end

        % 改，矩阵写法
        % 顺序sub_angle、s，分为4个部分
        for e_id_part = 1:4
            switch e_id_part
                case 1
                    % 与来流平行
                    e_id = (sub_angle_cos == 0);
                    cp(e_id) = 0;
                case 2
                    % 迎风面
                    e_id = (sub_angle_cos >0);
                    cp(e_id) = d * sub_angle_cos(e_id) .^2;
                case 3
                    % 背风面: Prandtl-Meyer method
                    e_id = (sub_angle_cos <0);
                    del = -(sub_angle(e_id) - pi/2);
                    cp(e_id) = -(del.^2)*(ga/2 + 1/2).*((16./(del.^2*Ma.^2*(ga + 1).^2) + 1).^(1/2) - 1);
                case 4
                    % 几何瑕疵导致面积为0
                    % cp结果是nan，一般是由于背风面计算中角度过小导致的当作与来流平行处理
                    e_id = (s ==0)|(isnan(cp));
                    cp(e_id) = 0;
            end            
        end



        e_force = cp.*-n1.*s/1e6;
        force_all_xyz = sum(e_force); % 气动合力
        force_pos = zeros(3,3);
        for force_xyz = 1:3
            % 求气动合力位置
            temp1 = p0; temp1(:,force_xyz) = 0;
            force_pos(force_xyz,:) = sum(e_force(:,force_xyz).*temp1)/force_all_xyz(force_xyz);
        end
        aero_pressure_center = sum(e_force(:,3).*p0(:,1))/force_all_xyz(3);% 气动合力位置，全机的压力中心：飞机上的总空气动力的作用线与飞机纵轴的交点称为
        aero_pressure_center_ba(i,j) = aero_pressure_center/L_r/1e3;
        % 力矩
        force_moment = cross(force_pos(1,:) - mass_center,[force_all_xyz(1) 0 0]);
        force_moment = force_moment + cross(force_pos(2,:) - mass_center,[0 force_all_xyz(2) 0]);
        force_moment = force_moment + cross(force_pos(3,:) - mass_center,[0 0 force_all_xyz(3)]);
        
        co_force_axis(i,j) = force_all_xyz(1)/S_r; % 轴向力系数
        co_force_lateral (i,j) = force_all_xyz(2)/S_r; % 侧向力系数
        co_force_normal(i,j) = force_all_xyz(3)/S_r; % 法向力系数
        co_force_lift(i,j) = co_force_axis(i,j)*sind(alpha) + co_force_normal(i,j)*cosd(alpha);% 升力系数
        co_force_drag(i,j) = (-co_force_axis(i,j)*cosd(alpha) + co_force_normal(i,j)*sind(alpha))*cosd(beta) + co_force_lateral(i,j)*sind(beta);% 阻力系数
        co_force_moment_roll(i,j) = force_moment(1)/S_r/L_r; % 滚转系数
        co_force_moment_pitch(i,j) = force_moment(2)/S_r/L_r;% 俯仰系数
        co_force_moment_yaw(i,j) = force_moment(3)/S_r/L_r;% 偏航系数
    end
end


% 压力系数分布图
% figure
% trisurf(Connectivity,Points(:,1),Points(:,2),Points(:,3),cp,'EdgeColor','none')
% axis equal
end

function s = tri_area_mat(con,point)
% 矩阵化
%% 根据空间中的三个点，计算面积
tri_point = zeros(size(con,1),3,3);
for i = 1:3
    tri_point(:,:,i) = point(con(:,i),:);
end
tri_len = zeros(size(con,1),3);
for i = 1:3
    temp = [i i+1];
    temp(temp>3) = 1;
    tri_len(:,i) = sqrt(sum((tri_point(:,:,temp(1))-tri_point(:,:,temp(2))).^2,2));
end
tri_p = sum(tri_len,2)/2;
s = sqrt(tri_p.*(tri_p-tri_len(:,1)).*(tri_p-tri_len(:,2)).*(tri_p-tri_len(:,3)));

end

