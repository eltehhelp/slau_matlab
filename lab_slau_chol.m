function [x,ok]=lab_slau_chol(A,b)
   ok = true;
   [m,n]=size(A);
 x=zeros(m,1); 
if  ~all(eig(A)>0)
    ok = false;    
         fprintf('ERROR: Chol method is not applicable.\n\n');
           disp(x)
       
else
L=zeros(m,n);
y=zeros(m,1);

for i=1:1:m
    for  j=1:i-1
        sum1=0;
    for   k=1:1:j-1
  sum1 =sum1+L(i,k)*L(j,k);
    end
    L(i,j)=(1/L(j,j))*(A(i,j)-sum1);
    end
    sum2=0;
     for k=1:1:i-1
         sum2=sum2+L(i,k)*L(i,k);
     end
      L(i,i)=sqrt(A(i,i)-sum2);
end

      %Ax=b/Ly=b
      
T=L';


for i=1:1:n
    sum=0;
    for k=i-1:-1:1
        
  sum=sum+y(k)*L(i,k) ;
    end    
  y(i)= (b(i)-sum)/L(i,i);
    
end



for i=n:-1:1
    sum=0;
    for k=i+1:1:n
        
    sum=sum+x(k)*T(i,k)  ;
    end
    
   x(i)= (y(i)-sum)/T(i,i);
end

%if show
 %   fprintf('Решение для метода Холецкого: \n');
  %  disp(x);
%end
 
end

