function gradient = loss_grad_function(L,x,del,alpha,S,anc)
if all(size(del)~=size(x))
    error('输入参数维度错误')
end
if size(anc,1) == size(x,1)
    anc_full = anc;
else
    anc_full = sparse(size(x,1),1);
    anc_full(logical(sum(S,2))) = anc;
end
term1 = 2*((L'*L*x)- L'*del);
term2 = 2* alpha *(S'*S*x - S'*anc_full);
gradient = term1 + term2;
end