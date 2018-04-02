function [h2_sd,h2_cp,h2_ss] = heregress(y,X)
% X - genotype matrix (N*M), where N is sample size and M is the number of
% markers.

ys=zscore(y);
ly1=i_extract_triu(i_square_diff_sum(ys,'diff'));   % he_sd
ly2=i_extract_triu(ys'*ys);                         % he_cp % y(1)*y(2) == 0.25*((y(1)+y(2))^2-(y(1)-y(2))^2) % he_cp mean corrected
ly3=i_extract_triu(i_square_diff_sum(ys,'sum'));    % he_ss

    f=0.5*nansum(X)./sum(~isnan(X));
    GRM=zeros(n);
    for i=1:n-1
        for j=i+1:n
            GRM(i,j)=0.5*nanmean((X(i,:)-2*f).*(X(j,:)-2*f)./(f.*(1-f)));
        end
    end
    lX=i_extract_triu(GRM);
    lXf=[ones(size(lX)) lX];
    b1=lXf\ly1;
    b2=lXf\ly2;
    b3=lXf\ly3;
    
    h2_sd=-b1(2)/b1(1);
    h2_cp=b2(2);
    h2_ss=b3(2);
end


function [X]=i_square_diff_sum(y,type)
    n=length(y);
    X=zeros(n);
    for i=1:n-1
        for j=i+1:n
            switch type
                case 'diff'
                    X(i,j)=(y(i)-y(j))*(y(i)-y(j));
                case 'sum'
                    X(i,j)=(y(i)+y(j))*(y(i)+y(j));                    
            end
        end
    end
end

function [x]=i_extract_triu(X)
    n=size(X,1);
    x=[];
    for k=1:n-1
        x=[x;diag(X,k)];
    end
end
