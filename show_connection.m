function label = show_connection(input)
unlabeled = 1:size(input.ConnectivityList,1);
label = zeros(size(input.ConnectivityList,1),1);
branch = 1;
imp = 0;
pos = 1;
label(1) = pos;
while any(label==0)
    new_branch = [];
    for i = 1:length(branch)
        for j = 1:3
            temp_idx = input.ConnectivityList(branch(i),j);
            new_branch = [new_branch;find(any(input.ConnectivityList==temp_idx,2)&(label==0))];
        end
    end
    if isempty(new_branch)
        pos = pos+1;
        branch = find(label==0,1);
    else
        label(new_branch) = pos;
        branch = new_branch;
    end
    disp([pos sum(label~=0)/length(label)])
end

end