function x = gseid(a,b,p,m,delta)
% Gauss-Seidel������
%Input	-a n��n���� -b n��1���� -p ��ʼ����n��1����
%   		-m ѭ������ -delta ����
%Output -x n��1����,Ax = b
n = length(b);
x = zeros(n,1);
for k = 1:m
    for j = 1:n
        if j ~= 1 && j ~= n
            x(j) = (b(j)-a(j,1:j-1)*x(1:j-1)...
                -a(j,j+1:n)*p(j+1:n))/a(j,j);
        elseif j == 1
            x(1) = (b(1)-a(1,2:n)*p(2:n))/a(1,1);
        else
            x(n) = (b(n)-a(n,1:n-1)*x(1:n-1))/a(n,n);
        end
    end
    err = abs(norm(x-p));
    relerr = err/(norm(x)+eps);
    p = x;
    if (err < delta) || (relerr< delta)
        break;
    end
end