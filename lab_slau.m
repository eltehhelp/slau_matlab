clc;
clear variables;
A=[1,-0.2589,-0.3093; -0.2589,1,-0.2705; -0.3093,-0.2705,1];
b=[1;1;1];
[m,n]=size(A);
show='true';
tol=1e-6;
e=10^(-10);
x0=zeros(m,1);
kmax=1000;
N=5000;

disp('������ ������ ���������� ������� ����');
[x,ok]=lab_slau_gauss(A,b);
fprintf('������� ��� ������ ������:\n');
           disp(x);
[x,ok] = lab_slau_kramer(A,b,tol); 
fprintf('������� ��� ������ �������:\n');
            disp(x);
[x,ok] = lab_slau_gauss_jordan(A,b);
fprintf('������� ��� ������ ������ - �������: \n');
    disp(x);
[x,ok]= lab_slau_matrix(A,b,tol);
fprintf('������� ��� ������ �������� �������:\n');
            disp(x);
[x,ok]=lab_slau_chol(A,b);
fprintf('������� ��� ������ ���������: \n');
    disp(x);

disp('������������ ������ ���������� ������� ����');
[x,ok,k] = lab_slau_jacobi(A, b, tol, x0, kmax );
fprintf('������� ��� ������ �����:\n');
            disp(x);
[x,ok,k]= lab_slau_gauss_seidel(A, b, tol, x0, kmax);
 fprintf('������� ��� ������ ������ �������:\n');
    disp(x);

[x,ok,k]=lab_slau_relax(A, b, e, x0, kmax);
fprintf('������� ��� ������ ����������:\n');
    disp(x);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load( 'lab_slau_data.mat', 'D'); %%��������� ������ ������ �� �����
x0=0;
s=size(D); 
tak=zeros(1,N);
T=zeros(1,8);
X=cell(size(D));
K=zeros(1,3);
X1=[];
for i=1:1:s(1,1)
    A=D{i}{1};
    b=D{i}{2};
    tak=zeros(1,N);
    for j=1:N
        tic %%������ �������(��������� ������� �����)
       [X1(:,1), ok]=lab_slau_gauss(A,b);
       tak(j)=toc; %%toc �������� �� �������� ������� �� ��� ���� ������������� ������
      if ~ok %%���� ����� �� ��������, �� ������� �� �����
          break
      end
    end
      
       T(i,1)=mean(tak)*ok; %%����������� �������������� ��������
    %%
    tak=zeros(1,N);
    for j=1:N
        tic
       [X1(:,2), ok]=lab_slau_kramer(A,b,tol); 
       tak(j)=toc;
      if ~ok
          break
      end
     end
      T(i,2)=mean(tak)*ok;
    %%
    tak=zeros(1,N);
    for j=1:N
        tic
      [ X1(:,3), ok]= lab_slau_gauss_jordan(A,b);
       tak(j)=toc;
      if ~ok
          break
      end
     end
    
       T(i,3)=mean(tak)*ok;
    %%
    tak=zeros(1,N);
    for j=1:N
        tic
      [ X1(:,4),ok]=lab_slau_matrix(A,b,tol);
          tak(j)=toc;
      if ~ok
          break
      end
    end
       T(i,4)=mean(tak)*ok;
    %%
    tak=zeros(1,N);
     for j=1:N
        tic
      [ X1(:,5),ok]=lab_slau_chol(A,b);
       tak(j)=toc;
      if ~ok
          break
      end
      end
      T(i,5)=mean(tak)*ok;
     
     %%
     tak=zeros(1,N);
    for j=1:N
        tic
       [X1(:,6),ok, K(i,1)]=lab_slau_jacobi(A, b, e, x0, kmax);
          tak(j)=toc;
      if ~ok
          break
      end
      T(i,6)=mean(tak)*ok;
    end
    %%
    tak=zeros(1,N);
    for j=1:N
        tic
       [X1(:,7), ok, K(i,2)]=lab_slau_gauss_seidel(A, b, e, x0, kmax);
        tak(j)=toc;
      if ~ok
          break
      end
    end
       T(i,7)=mean(tak)*ok;
      %%
      tak=zeros(1,N);
      for j=1:N
        tic
       [X1(:,8), ok, K(i,3)]=lab_slau_relax(A, b, e, x0, kmax);
        tak(j)=toc;
      if ~ok
          break
      end
       end
    T(i,8)=mean(tak)*ok;
      [m,n]=size(A);
   
 X{i}=X1; %%������ ������ �(� i�� �������� ��������� ������� ������� ���������� ����� �������� ��� i����)
  

%l=size(A);
r{i}=sprintf('A=%dx%d',m,n); %%����������� ������ �
  
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 fprintf('A=%dx%d double\n',m,n);
disp(A);
fprintf('b=%dx%d double\n',m,1);
disp(b);
%for i=1:1:s(1,1)
%s{i}=sprintf('A=%dx%d',m,n);
%end
%l=size(C);
 %fprintf('C= %d X %d \n', l);
 %disp(C);
 tablichka=table(X{i}(:,1), X{i}(:,2), X{i}(:,3), X{i}(:,4), X{i}(:,5),...
     'VariableNames',{'Gauss','Kramer','Gauss_Jordan', 'Matrix', 'Chol'});%% �������� �������
 disp(tablichka)%% ����� �������
figure( 'color', 'r')

subplot(2,2,1)
Tt=T';


bar(Tt(1:5, :))
grid on
grid minor
set(gca,'XTickLabel',{'Gauss','Kramer','Gauss Jordan','Matrix','Chol'});
ax=gca;
ax.Title.String={'������ ������','����� ����������'};
ax.Title.FontSize = 18;
ax.Title.Color = 'y';
ax.Title.FontName ='courier';
%title({'������ ������','����� ����������'});
legend(r,'location','northeastoutside');
subplot(2,2,2)



%c=categorical({'Jacobi', 'Gauss Seidel', 'Relax'});
bar(Tt(6:8, :))
grid on
grid minor
ax=gca;
set(gca,'XTickLabel',{'Jacobi','Gauss Seidel','Relax'});
ax.Title.String={'������������ ������','����� ����������'};
ax.Title.FontSize = 18;
ax.Title.Color = 'y';
ax.Title.FontName ='courier';
%title({'������������ ������','����� ����������'});
legend(r,'location','northeastoutside');

subplot(2,2,3)


H=Tt;
H(4,:)=[];

%c=categorical({'Gauss', 'Kramer', 'Gauss Jordan','Chol'});
bar(H(1:4,:))
grid on
grid minor
set(gca,'XTickLabel',{'Gauss','Kramer','Gauss Jordan','Chol'});
ax=gca;
ax.Title.String={'������ ������(��� ����������)','����� ����������'};
ax.Title.FontSize = 18;
ax.Title.Color = 'y';
ax.Title.FontName ='courier';
%title({'������ ������(��� ����������)','����� ����������'});
legend(r,'location','northeastoutside');

subplot(2,2,4)
%c=categorical({'Jacobi', 'Gauss Seidel', 'Relax'});
bar(K')
grid on
grid minor
ax=gca;
ax.Title.String={'������������ ������','���������� �������� '};
ax.Title.FontSize = 18;
ax.Title.Color = 'y';
ax.Title.FontName ='courier';
%title({'������������ ������','���������� �������� '})
set(gca,'XTickLabel',{'Jacobi','Gauss Seidel','Relax'});
legend(r,'location','northeastoutside');

   










