% graficos de salida del programa C++ de c´alculo de curvas principales
% Version que admite curvas recursivas
% Use:
%     grafsPCrec(datafile,outputfile);
%
% (C) Pedro Delicado (http://www-eio.upc.es/~delicado/PCOP)
% (C) Mario Huerta (http://www.revolutionresearch.org)
% 
function []=grafsPC(datafile,outputfile);

if (nargin==2)
   eval(['load ',datafile,'.dat']);
   eval(['X= ',datafile,';']);
   eval(['load ',outputfile,'.dat']);
   eval(['O= ',outputfile,';']);
else
   eval(['load ',datafile,'.dat']);
   eval(['X= ',datafile,';']);
   eval(['load ',datafile,'.alpha']);
   eval(['O= ',datafile,';']);
end    

prof=O(:,1);
O(:,1)=[];
% quitamos la coordenada de prof, y asi podemos usar el c´odigo de grafsPC.m
p0=(prof==0);
p1=(prof==1);

I=O(p0,1);
T=length(I);
p=length(O(1,:));
dim=(p-4)/2;

nc=100; % num. de puntos en los que se calcula el spline interpolador para los graficos
Ic=I(1)+(I(T)-I(1))*[0:(nc-1)]'/(nc-1);

density=O(p0,2);
span=O(p0,3);
var_k=O(p0,4);
alpha=O(p0,5:(5+dim-1));
b_ast=O(p0,(5+dim):p);

a_spl=zeros(nc,length(alpha(:,1)));

alpha1=O(p1,5:(5+dim-1));

%sum(density.*I.*diff([I(1);I])); % Esto es aproximadamente la media de los I, que debe ser aprox. 0

%graficos
figure(1)
lb=['density'; 'span   '; 'Var.res'];  
for i=1:3
   subplot(3,1,i)
   spl=spline(I,O(p0,(i+1)),Ic);  
   plot(Ic,spl,'-',I,O(p0,(i+1)),'.');
   title(lb(i,:))
end

figure(2)
for i=1:dim
   subplot(dim,1,i);
   a_spl(:,i)=spline(I,alpha(:,i),Ic);  
   plot(Ic,a_spl(:,i),'-',I,alpha(:,i),'.');
   title(['X_',num2str(i)])
end

switch dim
case 2
    figure(3)
    hold off
    plot(alpha(:,1),alpha(:,2),'ro')
    hold on
%   bfplot(alpha(:,1),alpha(:,2),'m',2)
    bfplot(a_spl(:,1),a_spl(:,2),'m',2)
%    for i=1:T
%       text(alpha(i,1),alpha(i,2),num2str(i))
%    end
    plot(X(:,1),X(:,2),'b.')
    xlabel('X_1')
    ylabel('X_2')
% case 3    
otherwise
    figure(3)
    hold off
    plot3(X(:,1),X(:,2),X(:,3),'b.')
    hold on
%   bfplot3(alpha(:,1),alpha(:,2),alpha(:,3),'m',2)
    plot3(alpha(:,1),alpha(:,2),alpha(:,3),'ro')
    plot3(alpha1(:,1),alpha1(:,2),alpha1(:,3),'go')
    plot3(alpha1(:,1),alpha1(:,2),alpha1(:,3),'b')
    bfplot3(a_spl(:,1),a_spl(:,2),a_spl(:,3),'m',2)
    xlabel('X_1')
    ylabel('X_2')
    zlabel('X_3')
end    

switch dim
case 2
    figure(4)
    hold off
    plot(alpha(:,1),alpha(:,2),'ro')
    hold on
    for i=1:T
       text(alpha(i,1),alpha(i,2),num2str(i))
    end
    plot(X(:,1),X(:,2),'b.')
    ejes=axis;
    eps=min(max(ejes)-min(ejes))/10; 
    for i=1:T
        plot([alpha(i,1)-eps*b_ast(i,1),alpha(i,1)+eps*b_ast(i,1)],...
             [alpha(i,2)-eps*b_ast(i,2),alpha(i,2)+eps*b_ast(i,2)],'r')
        plot(alpha(i,1)-eps*b_ast(i,1),alpha(i,2)-eps*b_ast(i,2),'+r')
    end
    xlabel('X_1')
    ylabel('X_2')
%case 3    
otherwise
    figure(4)
    hold off
    plot3(X(:,1),X(:,2),X(:,3),'b.')
    hold on
    plot3(alpha(:,1),alpha(:,2),alpha(:,3),'ro')
    plot3(alpha1(:,1),alpha1(:,2),alpha1(:,3),'go')
    plot3(alpha1(:,1),alpha1(:,2),alpha1(:,3),'b')
    ejes=axis;
    eps=min(max(ejes)-min(ejes))/10; 
    for i=1:T
        plot3([alpha(i,1)-eps*b_ast(i,1),alpha(i,1)+eps*b_ast(i,1)],...
             [alpha(i,2)-eps*b_ast(i,2),alpha(i,2)+eps*b_ast(i,2)],...
             [alpha(i,3)-eps*b_ast(i,3),alpha(i,3)+eps*b_ast(i,3)],'r')
        plot3(alpha(i,1)-eps*b_ast(i,1),alpha(i,2)-eps*b_ast(i,2),alpha(i,3)-eps*b_ast(i,3),'+r')
    end
    xlabel('X_1')
    ylabel('X_2')
    zlabel('X_3')
end    

