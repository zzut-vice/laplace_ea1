% 超松弛迭代法
function x = sor(A,b,x0,w,tol)
D=diag(diag(A));
L = -tril(A,-1);
U = -triu(A,1);
B = (D - w * L) \ ((1 - w) * D + w * U);
g = (D - w * L) \ b * w;
x = B * x0 +g;
er = norm(x - x0);
n=1;
while er >= tol
    x0 = x;
    x= B * x0 + g;
    n = n + 1;
    er = norm(x - x0);
    disp(er)
end
end
