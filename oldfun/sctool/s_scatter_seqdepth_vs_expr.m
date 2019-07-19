load examgrades.mat

%%
[~,i]=sort(sum(X,2),'descend');
X=X(i,:);

figure;
hold on
%subplot(2,2,1)
i=1:2500;
scatter(sum(X(i,:)),sum(X(i,:)>0));
refline;

%subplot(2,2,2)
i=2501:5000;
scatter(sum(X(i,:)),sum(X(i,:)>0));
refline;

%subplot(2,2,3)
i=5001:7500;
scatter(sum(X(i,:)),sum(X(i,:)>0));
refline;

%subplot(2,2,4)
%scatter(sum(X(7500:end,:)),sum(X(7500:end,:)>0));
refline;




%%

figure;
hold on
%subplot(2,2,1)
i=1:2500;
scatter(sum(X),sum(X(i,:)>0));
refline;

%subplot(2,2,2)
i=2501:5000;
scatter(sum(X),sum(X(i,:)>0));
refline;

%subplot(2,2,3)
i=5001:7500;
scatter(sum(X),sum(X(i,:)>0));
refline;

%subplot(2,2,4)
scatter(sum(X),sum(X(7500:end,:)>0));
refline;


%%
figure; hold on;
for k=1:100
scatter(log(sum(X)),log(X(k,:)),'k')
end


