function [x, ok]=lab_slau_gauss(A,b)
[m, n] = size(A);
   ok = true;
    x = zeros(m, 1);
if m~=n
    ok = false;
        
           fprintf('ERROR: Gauss method is not applicable.\n\n');
           disp(x)
        
else

for i=1:1:(m-1)
    for k=i+1:1:m%!! МАСЛЕННИКОВ СКАЗАЛ ТУТ БУДЕТ КОСЯК
        M=A(k,i)/A(i,i);
        A(k,:)=A(k,:)-A(i,:)*M ;    
b(k)=b(k)-b(i)*M;
 
    end
end

for i=n:-1:1
    sum=0;
    for k=i+1:1:n
        
    sum=sum+x(k)*A(i,k)  ;
    end
    
   x(i)= (b(i)-sum)/A(i,i);
   
    
end
        %if show
         %   fprintf('Решение для метода Гаусса:\n');
          %  disp(x);
       % end

end

