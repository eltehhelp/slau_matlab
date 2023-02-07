function [x,ok]=lab_slau_matrix(A,b,tol)
[m, n] = size(A);
   ok = true;
   % x = zeros(m, 1);
if abs(det(A))<tol
    ok = false;
        
           fprintf('ERROR: Kramers method is not applicable.\n\n');
           disp(x)
        
else
for i=1:m
    for j=1:n
 M=A;
 M(i,:)=[];
 M(:,j)=[];
 L=det(M);
%  k=;
B(i,j)=(-1)^(i+j).*L;
    end
end
B=B';
R=(1/det(A))*B;
%x=R*b
%A*x=b
x=R*b;
end

       % if show
        %    fprintf('Решение для метода обратной матрицы:\n');
         %   disp(x);
        %end
end

