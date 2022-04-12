function x=GaussianSolver(A,b)
[n,~] = size(A);
x=zeros(n,1);
 
for j = 1:n-1
    for i = j+1:n
        mul = A(i,j)/A(j,j);
        A(i,:)= A(i,:)-mul*A(j,:);
        b(i)= b(i)-mul*b(j);
    end
end
 
 
for i=n:-1:1
    sum=0;
    for j=n:-1:i+1
        sum=sum+x(j)*A(i,j);
    end
    x(i)=(b(i)-sum)/A(i,i);
end
end