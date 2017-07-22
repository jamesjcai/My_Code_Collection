function [Result,Value,Time]=epiACO(FileName,AntNumber,IterationNumber,evaporation)
%% EpiACO main program for the detection of SNP-SNP interactions
% Input
%   FileName:        The input file name. for example,
%                    FileName='dat100' or FileName='Model1_Set1'; Be
%                    sure that the file format should like the example file.
%   AntNumber:       The ant number. for example, AntNumber=50.
%   IterationNumber: The iteration number. for example, IterationNumber=20.
%   evaporation:     The evaporation of Tau. for example,evaporation=0.9.
%
% Output
%   Result:         The final results.
%   Value:          The Svalues of final results.
%   Time:           The computational time.
%
% Example
%                   [Result,Value,Time]=EpiACO('Model1_Set1',250,20,0.9)

%% Time Begin
tic;

%% Unimportant Parameters
% The considered order of SNP-SNP interactions.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ratio
ratio=0.5;

% The extration of Tau.
extration=0.3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Dimension=2;

% The initialized Tau value.
Tau0=1;

% The coefficient of Tau.
Alfa=1;

% The initialized Eta value.
Eta0=1;

% The coefficient of Eta.
Beta=1;

% The amplification of Svalue.
Amplification=10000;

%% Input data
tmp=load([FileName,'.mat']);
pts=tmp.pts;
class=tmp.class;
SNPNumber=size(pts,2);

%% Initialize Parameters
% Initialize Tau
Tau=Tau0*ones(1,SNPNumber);

% Initialize Eta
Eta=Eta0*ones(1,SNPNumber);

% Initialize Combination
Combination=zeros(AntNumber,Dimension);

% Initialize FitnessValue
FitnessValue=zeros(AntNumber,1);

% initialize Result
Result=[];

% initialize Value
Value=[];

%% Iteration

disp('Waitting...');

for i=1:IterationNumber
    %%
    % Select combinations and evaluate them.
    q0=i/IterationNumber;
    for j=1:AntNumber
        Combination(j,:)=SelectCombination(Tau,Eta,Alfa,Beta,q0,Dimension,ratio);
        FitnessValue(j,1)=Svalue(pts,class,Combination(j,:),Amplification);
    end
    
    %%
    % Rank combinations and Select suspected combinations.
    TempCombination=[Combination;Result];
    TempFitnessValue=[FitnessValue;Value];
    
    [TempCombination,RowIndex,~]=unique(TempCombination,'rows');
    TempFitnessValue=TempFitnessValue(RowIndex,:);
    
    [TempFitnessValue,Index]=sort(TempFitnessValue,'descend');
    TempCombination=TempCombination(Index,:);
    
    Tag=InflexionPoint(TempFitnessValue);
    Result=TempCombination(1:Tag,:);
    Value=TempFitnessValue(1:Tag,:);
    
    %%
    % Update pheromone
    Tau=UpdatePheromone(Tau,Combination,FitnessValue,Result,Value,evaporation,extration);
    
    %%%%%%%%%%%%%%%%%%
    if sum(Tau>mean(Tau))<=max(2,floor(SNPNumber*0.1))
        break;
    end
    %%%%%%%%%%%%%%%%%%
    
%     subplot(2,1,1);
%     plot(TempFitnessValue);
%     title('Svalue');
%     subplot(2,1,2);
%     plot(1:SNPNumber,Tau,'+b-.') ;
%     hold on
%     plot(1:SNPNumber,mean(Tau),'pr')
%     title('Tau');
%     hold off
%     drawnow;
    
end
% post processing
[Tau,Index]=sort(Tau,'descend');
Lab=InflexionPoint(Tau');

%%%%%%%%%%%%%%%%%%%%%%%%%%%
Lab=min([Lab,sum(Tau>mean(Tau)),5]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%

matrix=Index(1,1:Lab);
TempCombination1=nchoosek(matrix,Dimension);
TempCombination1=sort(TempCombination1,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%
TempCombination1=setdiff(TempCombination1,Result,'rows');
%%%%%%%%%%%%%%%%%%%%%%%%%%%
r=size(TempCombination1,1);
TempFitnessValue1=zeros(r,1);
for i=1:r
    TempFitnessValue1(i,1)=Svalue(pts,class,TempCombination1(i,:),Amplification);
end

TempCombination2=[Result;TempCombination1];
TempFitnessValue2=[Value;TempFitnessValue1];

[TempCombination2,RowIndex,~]=unique(TempCombination2,'rows');
TempFitnessValue2=TempFitnessValue2(RowIndex,:);

[TempFitnessValue2,Index]=sort(TempFitnessValue2,'descend');
TempCombination2=TempCombination2(Index,:);

Tag=InflexionPoint(TempFitnessValue2);
if Tag>5
    Tag=5;
end
Result=TempCombination2(1:Tag,:);
Value=TempFitnessValue2(1:Tag,:);

%% Time End
Time=toc;

%% Save Results
SaveName=[FileName,'_',num2str(AntNumber),'_',num2str(IterationNumber),...
    '_',num2str(evaporation),'.mat'];
save(SaveName,'Result','Value','Time');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Combination=SelectCombination(Tau,Eta,Alfa,Beta,q0,Dimension,ratio)
%% Select a SNP combination according to the probability formula
[~,index]=sort(Tau);
index1=index(1,1:floor(ratio*size(index,2)));
Number=size(index1,2);
Combination=zeros(1,Dimension);
if rand>q0
    combine=randperm(Number,Dimension);
    Combination=index(combine);
else
    for i=1:Dimension
        Weight=(Alfa*Tau).*(Beta*Eta);
        SumWeight=sum(Weight);
        Ratio=Weight/SumWeight;
        Combination(1,i)=Roulette(Ratio);
        Tau(1,Combination(1,i))=0;
    end
end
Combination=sort(Combination);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SNP=Roulette(Ratio)
%% A roulette

Sel=rand;
SumPs=0;
i=1;
while SumPs<Sel
    SumPs=SumPs+Ratio(1,i);
    i=i+1;
end
SNP=i-1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function FitnessValue=Svalue(pts,class,factor,Amplification)
%% fitness function Svalue

FitnessValue=Amplification*MI(pts,class,factor)*(1/BN(pts,class,factor));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function FitnessValue=MI(pts,class,factor)
%% Mutual Information value

data=pts(:,factor);
FitnessValue=Func_MutualInfo(double(data)',double(class));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function FitnessValue=Func_MutualInfo(data,labels)
% Jan.19 2007 By Guoqiang Yu
%
%--- INPUT
% data:  features x samples; note that the
% value of zero indicates a missing data
% labels: the label matrix of each sample;
%
%--- OUTPUT
% MI: Information provided by the features
%%
data=data-min(data(:));          %minimum data is 0
labels=labels-min(labels(:))+1; %minimum label is 1

Num_Label=max(labels(:));
Num_DataType=max(data(:))+1;
[Num_Feature,Num_Sample]=size(data);

%Entropy of the random variable Label
H_Label=hist(labels,Num_Label);
P_Label=H_Label/Num_Sample;
Entropy_Label=-P_Label*log2(P_Label'+eps);

%Special dealing with the case of small Num_Feature
if Num_Feature<9
    ZZ=Num_DataType.^(Num_Feature-1:-1:0)';
    Hist_Label=zeros(Num_DataType^Num_Feature,Num_Label);%  Hist_Label is p(c,f)
    tempIndex=ZZ'*data+1;
    for j=1:Num_Sample
        Hist_Label(tempIndex(j),labels(j))=Hist_Label(tempIndex(j),labels(j))+1; % calculate p(c,f)
    end
    
    sumHist=sum(Hist_Label,2);   %calculate p(f)
    repHist=repmat(sumHist,1,Num_Label);
    pHist_Label=Hist_Label./(repHist+eps);%p(c/f)=p(c,f)/p(f).
    InfoIncre=-sum((log2(pHist_Label+eps).*pHist_Label).*(repHist));
    FitnessValue=Entropy_Label-sum(InfoIncre)/Num_Sample;
    return;
end

%Larger Feature Number follows the following procedure
mm=1;
Hist_Label=zeros(Num_Label,Num_Sample);
Hist_Label(labels(1,1),mm)=1;
Hist_SNP=zeros(Num_Feature,Num_Sample);
Hist_SNP(:,mm)=data(:,1);

for j=2:Num_Sample
    tempData=data(:,j);
    Index=0;
    for k=1:mm
        if isequal(Hist_SNP(:,k),tempData)
            Index=k;
            break;
        end
    end
    if Index==0
        mm=mm+1;
        Hist_SNP(:,mm)=tempData;
        Hist_Label(labels(j,1),mm)=1;
    else
        Hist_Label(labels(j,1),Index)=Hist_Label(labels(j,1),Index)+1;
    end
end

M1=mm;
InfoIncre=0;
for s=1:M1
    tempNum=sum(Hist_Label(:,s));
    P=Hist_Label(:,s)/tempNum;
    InfoIncre=InfoIncre-P'*log2(P+eps)*tempNum;
end
FitnessValue=Entropy_Label-InfoIncre/Num_Sample;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function FitnessValue=BN(pts,class,factor)
%% The Bayes Network value

snp_com=pts(:,factor)-1;
state=class'-1;
[xrow,~] = size(snp_com);
subs = snp_com+1;
sample = accumarray(subs,ones(xrow,1));
disease = accumarray(subs,state);
control = sample-disease;
sample(4,4) = 0;
disease(4,4) = 0;
control(4,4) = 0;
z=0;
for i = 1:3
    for j = 1:3
        y=My_factorial(sample(i,j)+1);
        r=My_factorial(disease(i,j))+My_factorial(control(i,j));
        z=z+(r-y);
    end
end
FitnessValue=abs(z);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f=My_factorial(e)
f=0;
if e>0
    for o=1:e
        f=f+log(o);
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Tag=InflexionPoint(TempFitnessValue)
%% Find a inflexion point
PointNumber=size(TempFitnessValue,1);

Slope=zeros(1,PointNumber-1);
Slope2=zeros(1,PointNumber-2);

for i=1:(PointNumber-1)
    Slope(1,i)=TempFitnessValue(i,1)-TempFitnessValue(i+1,1);
end

for i=1:(PointNumber-2)
    Slope2(1,i)=abs(Slope(1,i)-Slope(1,i+1));
end

[~,Tag]=max(Slope2);
Tag=Tag+2;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Tau=UpdatePheromone(Tau,Combination,FitnessValue,Result,Value,evaporation,extration)
%% Update pheromone

DeltaTau=zeros(1,size(Tau,2));
for i=1:size(Combination,1)
    DeltaTau(1,Combination(i,:))=DeltaTau(1,Combination(i,:))+FitnessValue(i,1);
end
DeltaTau1=zeros(1,size(Tau,2));
for i=1:size(Result,1)
    DeltaTau1(1,Result(i,:))=DeltaTau1(1,Result(i,:))+extration*Value(i,1);
end
Tau=Tau.*(1-evaporation)+DeltaTau+DeltaTau1;
end
