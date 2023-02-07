function [x,ok] = lab_slau_gauss_jordan(A,b)
ok=true;
[m, n]=size(A);
if m ~= n 
    ok = false;
    
        fprintf('ERROR: GAUS_JORDAN is not applicable \n\n')
        disp(x)
    
else

C=[A,b];

for i=1:1:(m-1)
   for k=(i+1):1:m
   C(k,:)=C(k,:)-C(i,:)*C(k,i)/C(i,i);
   end
end


for i=m:(-1):1
for k=(i-1):-1:1
     C(k,:)=C(k,:)-C(i,:)*C(k,i)/C(i,i);
    
end
end

for i=1:m
   x(i)=(C(i,n+1)/C(i,i));
end
x=x';
%if show
 %   fprintf('Решение для метода Гаусса - Жордана: \n');
  %  disp(x);
%end
end


   

