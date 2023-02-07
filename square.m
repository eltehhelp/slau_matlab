%% метод квадратного корня для решения СЛАУ
function [a,b] = square(a,b)
n=length(b);
u(1,1)=sqrt(a(1,1));
for j=2:n
    u(1,j)=a(1,j)/u(1,1);
end
for i=2:n
    u(i,i)=sqrt(a(i,i)-sum(u(1:i-1,i).*u(1:i-1,i)));
    for j=i+1:n
        u(i,j)=(a(i,j)-sum(u(1:i-1,i).*u(1:i-1,j)))/u(i,i);
    end
end
y= u'\b;
x= u\y
end