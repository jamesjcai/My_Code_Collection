% Set Program Parameters
outtag='PANDAOutput'; % used to generate the output file-name
alpha=0.1; % update parameter, recommend values between 0.05 and 0.25
exp_file='YeastData/YeastData_KOExpression.txt'; % text file with expression values
motif_file='YeastData/YeastData_Motif.txt'; % text file with initial edges and weights for regulatory network
ppi_file=''; % text file with initial edges and weights for the TF-interaction network. Note that setting ppi_file='' initializes the PPI to an identity matrix

%% Read in Data %%
disp('Reading in data!')

% Expression Data
fid=fopen(exp_file, 'r');
headings=fgetl(fid);
while(strcmp(headings(1), '#'))
	headings=fgetl(fid);
end
NumConditions=length(regexp(headings, '\t'));
frewind(fid);
Exp=textscan(fid, ['%s', repmat('%f', 1, NumConditions)], 'delimiter', '\t', 'CommentStyle', '#');
fclose(fid);
GeneNames=Exp{1};
NumGenes=length(GeneNames);
Exp=cat(2, Exp{2:end});

% Prior Regulatory Network
[TF, gene, weight]=textread(motif_file, '%s%s%f');
TFNames=unique(TF);
NumTFs=length(TFNames);
[~,i]=ismember(TF, TFNames);
[~,j]=ismember(gene, GeneNames);
RegNet=zeros(NumTFs, NumGenes);
RegNet(sub2ind([NumTFs, NumGenes], i, j))=weight;

% PPI Data
TFCoop=eye(NumTFs);
if(~isempty(ppi_file))
	[TF1,TF2,weight]=textread(ppi_file, '%s%s%f');
	[~,i]=ismember(TF1, TFNames);
	[~,j]=ismember(TF2, TFNames);
	TFCoop(sub2ind([NumTFs, NumTFs], i, j))=weight;
	TFCoop(sub2ind([NumTFs, NumTFs], j, i))=weight;
end

% Build Co-expression Matrix
NumConditions=size(Exp,2);
if(NumConditions>2)
	GeneCoReg=corr(Exp', 'type', 'pearson', 'rows', 'pairwise');
	GeneCoReg(1:NumGenes+1:NumGenes^2)=1;
	GeneCoReg(isnan(GeneCoReg))=0;
else
	GeneCoReg=eye(NumGenes);
end

%% Run Message-Passing
AgNet=PANDA(RegNet, GeneCoReg, TFCoop, alpha);

% reshape network information into vectors and print to file
TF=repmat(TFNames, 1, length(GeneNames));
gene=repmat(GeneNames', length(TFNames), 1);
TF=TF(:);
gene=gene(:);
RegNet=RegNet(:);
AgNet=AgNet(:);

fid=fopen([outtag, '_FinalNetwork.pairs'], 'wt');
fprintf(fid, 'TF\tgene\tMotif\tPANDA-prediction\n');
for(cnt=1:length(TF))
	fprintf(fid, '%s\t%s\t%f\t%f\n', TF{cnt}, gene{cnt}, RegNet(cnt), AgNet(cnt));
end
fclose(fid);
