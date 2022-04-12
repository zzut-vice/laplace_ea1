function V = cal_volume(Connectivity,Points)
V = 0;
for i=1:length(Connectivity)
    points_temp = Points(Connectivity(i,:),:);
    vi=(-points_temp(3,1)*points_temp(2,2)*points_temp(1,3)...
        +points_temp(2,1)*points_temp(3,2)*points_temp(1,3)...
        +points_temp(3,1)*points_temp(1,2)*points_temp(2,3)...
        -points_temp(1,1)*points_temp(3,2)*points_temp(2,3)...
        -points_temp(2,1)*points_temp(1,2)*points_temp(3,3)...
        +points_temp(1,1)*points_temp(2,2)*points_temp(3,3))/6;
    V=V+vi;
    
%     xi1=vout(i,1);     yi1=vout(i,2);     zi1=vout(i,3);
%     xi2=vout(i+1,1);   yi2=vout(i+1,2);   zi2=vout(i+1,3);
%     xi3=vout(i+2,1);   yi3=vout(i+2,2);   zi3=vout(i+2,3);
%     vi=(-xi3*yi2*zi1+xi2*yi3*zi1+xi3*yi1*zi2-xi1*yi3*zi2-xi2*yi1*zi3+xi1*yi2*zi3)/6;
%     V=V+vi;
end
V = abs(V);

end