function aero_obj = aerofun(Points,ConnectivityList)
[co_force_lift,co_force_drag]...
    = CalPressure_alphama(Points,ConnectivityList,[10 15 20],0,[6 8 10]);
temp1 = co_force_lift./co_force_drag;
% 基准值
baseline = [1.75331784013013	1.65436590287488	1.59176582260317
1.89957371216340	1.85215500486033	1.82312234834351
1.81211614616233	1.79453741110935	1.78398294475615];
temp1 = (temp1 - baseline)./baseline;
aero_obj =  min(temp1(:));% 此处目标是各项中最差的提升百分比

end