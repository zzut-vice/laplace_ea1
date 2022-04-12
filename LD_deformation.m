function Points_df = LD_deformation(moving_points,...
    L_x1_inv,L_x2_inv,L_x3_inv,...
    L_x1,L_x2,L_x3,...
    del_x1_s,del_x2_s,del_x3_s)
% L_x'*L_x*x = L_x'*del_x
del_x1 = [del_x1_s;moving_points(:,1)];
Points_1 = L_x1_inv*L_x1'*del_x1;
del_x2 = [del_x2_s;moving_points(:,2)];
Points_2 = L_x2_inv*L_x2'*del_x2;
del_x3 = [del_x3_s;moving_points(:,3)];
Points_3 = L_x3_inv*L_x3'*del_x3;
Points_df = [Points_1 Points_2 Points_3];
end