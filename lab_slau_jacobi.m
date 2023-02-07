function [x,ok,k]=lab_slau_jacobi(A, b, e, x0, kmax)
[m,n]=size(A);
D=diag(diag(A));
    B=(D^(-1))*(D-A);
    ok=true;
    q=norm(B);
   k=0;
   x=zeros(m,1);
   
    if  ~all(abs(eig(B))<1)
        ok = false;
        
            fprintf('ERROR: Jakobi method is not applicable.\n\n');
            disp(x);
        
                
    else
    x=ones(m,1);
      

    while (k<=kmax)&&~(max(abs(x-x0))/(1-q)<=e)
        k=k+1;
        x0=x;
    for i=1:1:m
        sum=0;
        for j=1:1:n
            if(j~=i)
               sum=sum+A(i,j)*x0(j); 
            end
        end
        x(i)=(b(i)-sum)/A(i,i);
    end
    
    end
 
    end
end


