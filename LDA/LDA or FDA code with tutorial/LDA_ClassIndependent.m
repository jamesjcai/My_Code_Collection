% Prophet Mohammed said [ALLAH will help any one helped his/her brother/sister] PBUH

function [W,y] = LDA_ClassIndependent(Data, a)

% This code is written by Alaa Tharwat (Frankfurt University of Applied Science - Germany)
% For more details about the code of the numerical example(s) see our paper "Tharwat, A., Gaber, T., Ibrahim, A., 
% & Hassanien, A. E. Linear discriminant analysis: A detailed tutorial. AI Communications,
% (Preprint), 1-22.?
% engalaatharwat@hotmail.com

%% Examples
% c1=[1 2;2 3;3 3;4 5;5 5]  % the first class (5 observations)
% c2=[4 2;5 0;5 2;3 2;5 3;6 3] % the second class (6 observations)
% data=[c1;c2] % te whole data
% a=[5;6] % the number of samples in each class
% [W,y] = LDA_ClassIndependent(data, a)


% Calculate the total number of all classes
c=unique(a);
pos=1;
for i=1:length(c)
   for j=1:a(i)
      Labels(pos,1)=i;
      pos=pos+1;
   end
end

% Calculate Mean of each class
for i=1:length(c)
    mu(i,1:size(Data,2))=mean(Data(Labels==i,:));
end
% Calculate the total mean of all classes
muTotal=zeros(size(mu(1,:)));
for i=1:length(c)
    muTotal=muTotal+a(i)*mu(i,:);
end
muTotal=muTotal/(sum(a));

% Subtract the originla data from the mean
D=zeros(size(Data,1),size(Data,2));
for i=1:length(c)
    D(Labels==i,:)=Data(Labels==i,:)-repmat(mu(i,:),a(i),1);
end

% Calculate the within class variance (SW)
SW=zeros(size(Data,2),size(Data,2));
for i=1:c
    SW=SW+D(Labels==i,:)'*D(Labels==i,:);
end

% Calculate the Between-class variance (SB)
SB=zeros(size(Data,2),size(Data,2));
for i=1:length(c)
   SB= SB+a(i)*(mu(i,:)-muTotal)'*(mu(i,:)-muTotal);
end

% Calculate J(W)
J=inv(SW)*SB;

% Calculate the eignevalues and eigenvectors of (J)
[evec,eval]=eig(J);

% Sort the eigenvectors according to their corresponding eigenvalues (descending order)
eval = diag(eval);
[junk, index] = sort(-eval);
eval = eval(index);
evec = evec(:, index);

% Slect the most largest c eigenvectors as a lower dimensional space
W=evec(:,1:length(c)-1);
% project the original data on thelower dimensional space (W)
y=Data*W;
