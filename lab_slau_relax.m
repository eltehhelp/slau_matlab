function [x, ok,k] =lab_slau_relax(A, b, e, x0, kmax, show)
[m,n]=size(A);
L=tril(A,-1);
U=triu(A,1);
D=diag(diag(A));
P=-(L+D)^(-1)*U;
w=1.1;
ok=true;

    if ~all(abs(eig(P))<1)
        ok = false;
         
            fprintf('ERROR: relax method is not applicable.\n\n');
            disp(x)
        
       
    else
x=ones(m,1);
    k=0;
  while (k<=kmax) && (norm(A*x0-b)>=e) 
      x0=x;
    for i=1:1:m
        sum1=0;
        sum2=0;
        for j=1:1:i-1
            sum1=sum1+A(i,j)*x(j); 
        end
        for j=i+1:1:n
            sum2=sum2+A(i,j)*x0(j); 
        end
       x(i)=(1-w)*x0(i)+w/A(i,i)*(b(i)-sum1-sum2);
    end
     k=k+1;
  end
   % if show
    %fprintf('Решение для метода релаксации:\n');
    %disp(x);
    %end
    end
   
end
