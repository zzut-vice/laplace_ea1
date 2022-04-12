%% Gauss-Serdel迭代法的函数文件gauseidel.m：  
function [y,n]=gauseidel(A,b,x0,ep)  
D=diag(diag(A));  
L=-tril(A,-1);  
U=-triu(A,1); 
B=(D-L)\U;  
f=(D-L)\b;  
y=B*x0+f;  
er = norm(y-x0);
n=1;  
while er>=ep  
    x0=y;
    y=B*x0+f;
    n=n+1;
    er = norm(y-x0);
    disp(['dif:' num2str(er)])
end
end