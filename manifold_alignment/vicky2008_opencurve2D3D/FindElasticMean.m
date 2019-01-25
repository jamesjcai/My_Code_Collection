function [p,mu,q,E] = FindElasticMean(Data)

Niter= 20;
n=size(Data,3); 
figure(1);clf; 
for i=1:n
%     X = ReSampleCurve(Data(:,:,i),200);
    q(:,:,i) = curve_to_q(Data(:,:,i));
%     q(:,:,i) = ProjectC(q(:,:,i));
end

del = .5;
mu =  q(:,:,3);
E(1)=1;
iter=2;

while (iter<=Niter && E(iter-1)>0.05)
    vm = 0;    
    for i=1:n
        [iter i]
        v = ElasticShootingVector(q_to_curve(mu),Data(:,:,i),mu,q(:,:,i),1);
        vm = vm + v;

    end
       vm = vm/n;

       E(iter) = norm(vm,'fro');
       mu = ElasticShooting(mu,del*vm);
%        mu = ProjectC(mu);
       iter=iter+1;
       E   
end



figure(11);clf; 
plot(E(2:end));


figure(21); clf; 
p = q_to_curve(mu);
plot3(p(1,:),p(2,:),p(3,:),'LineWidth',3);axis equal;axis off;
% plot(p(1,:),p(2,:),'LineWidth',3)
axis equal;
% %a1=C(:,:,1);a2=C(:,:,5);a3=C(:,:,9);a3=C(:,:,13);a4=C(:,:,17);a4=C(:,:21);a5=C(:,:,25);a6=C(:,:,29);a7=C(:,:,33);a8=C(:,:,37);a9=C(:,:,41);a10=C(:,:,45);
% Data=(a1;a2;a3)
% Data=A
% for i=1:5:50 
% A(:,:,length(i))=C(:,:,i);
% end
% % figure(i);
% % subplot()
%clear;close all;A=C(:,:,1:5:50);Data=A(:,1:70,:);[p,mu,q,E] = FindElasticMean(Data);
