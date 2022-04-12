function loss = loss_function(L,x,del,alpha,S,anc)
if all(size(del)~=size(x))
    error('输入参数维度错误')
end
shape_term_vec = (L*x-del);
shape_term = sum(shape_term_vec.^2,'all');

if size(anc,1) == size(x,1)
    anc_full = anc;
else
    anc_full = sparse(size(x,1),1);
    anc_full(logical(sum(S,2))) = anc;
end

con_term_vec = alpha*(S*x-anc_full);
con_term = sum(con_term_vec.^2,'all');

loss = shape_term + con_term;


end