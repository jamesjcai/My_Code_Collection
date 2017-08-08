function i_joyplot3col(D,g,scalek,styleid)

% James Cai (jcai@tamu.edu)
% (c) Aug 2017


    if nargin<4, styleid=1; end
    if nargin<3 || isempty(scalek), scalek=0.3; end
    if nargin<2 || isempty(g), g=zeros(size(D,1),1); end

    D=D./max(D(:));
    a=max([sum(g==0), sum(g==1) sum(g==2)]);
    D=D(:,1:50:end);
    for kx=1:3
        switch styleid
            case 1
                hold on
                i_joy(D(g==kx-1,:),scalek,kx);
            otherwise
                hold on
                i_joy2(D(g==kx-1,:),scalek,kx);
        end
        ylim([0 (a*scalek)+(1-scalek)]);
        axis off
    end
end

function i_joy(D,scalek,coln)
    if nargin<3, coln=1; end
    [n,p]=size(D);
    for k=n:-1:1
       plot(1.1*p*(coln-1)+(1:p),...
            scalek*(k-1)*ones(1,p)+D(k,:),'k-');
    end
    xlim([1 3*p+0.2*p]);
    set(gca, 'YTick', []);
end

function i_joy2(D,scalek,coln)
    if nargin<3, coln=1; end
    [n,p]=size(D);
    x=1:(p+2);
    y=zeros(size(x));
    for k=n:-1:1
       y(1)=scalek*(k-1);
       y(2:end-1)=scalek*(k-1)*ones(1,p)+D(k,:);
       y(end)=scalek*(k-1);
       x1=1.1*(p+2)*(coln-1)+x;
       patch(x1,y,'red','FaceAlpha',.5);
    end
    xlim([1 3*p+0.2*p]);
    set(gca, 'YTick', []);
end

% function i_joy3(D)
% [n,p]=size(D);
%  for k=1:n
%      subplot_tight(n, 1, k, 0);    
%      bar(D(k,:),1);
%      xlim([1 p]);
%      ylim([0 1]);
%      set(gca, 'XTick', []);
%      set(gca, 'YTick', []);
%  end
% end

