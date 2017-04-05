function Matrix_eQTL_engine(snps,gene,cvrt,Output_file_name,pvOutputThreshold,useModel,errorCovariance, verbose)
% Matrix eQTL by Andrey A. Shabalin
% http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/
% Build 1.2.0

% Set missing parameters
if(nargin<8)
	verbose = false;
	if(nargin<7)
		errorCovariance = [];
		if(nargin<6)
			useModel = modelLINEAR;
			if(nargin<5)
				pvOutputThreshold = 1e-5;
				if(nargin<4)
					error('Too few parameters');
				end;				
			end;
		end;
	end;
end;

gene = gene.Clone();
%snps = snps.Clone();
cvrt = cvrt.Clone();

switch useModel
	case modelLINEAR
		statistic_name = 't-stat';
		snps_process = @impute_row_mean;
		nVarTested = 1;
	case modelANOVA
		snps_process = @snps_split_for_ANOVA;
		statistic_name = 'F-test';
		nVarTested = 2;
	otherwise
		error('Unknown value of useModel');
end;

if(verbose)
	tic
	status = @status_T;
else
	status = @status_F;
end;

% Check dimensions
status('Checking input data dimensions',[]);
if(snps.nCols()*snps.nRows() == 0)
	error('Empty genotype dataset');
end;
if(gene.nCols()*gene.nRows() == 0)
	error('Empty expression dataset');
end;
if(snps.nCols ~= gene.nCols)
	error('Different number of samples in the genotype and gene expression files');
end;
if(cvrt.nRows>0)
	if(snps.nCols ~= cvrt.nCols)
		error('Wrong number of samples in the file with covariates');
	end;	
end;
%% ################################# error covariance processing #################################

if(~isempty(errorCovariance))
	status('Processing the errorCovariance matrix');
	if(size(errorCovariance,1)~=size(errorCovariance,2))
		error('The covariance matrix is not square');
	end;
	if(size(errorCovariance,1)~=snps.nCols)
		error('The covariance matrix size does not match the data');
	end;	
	% test for symmetry
	if(~all(all(errorCovariance==errorCovariance')))
		error('The covariance matrix is not symmetric');
	end;
	[v,d] = eig(errorCovariance);	
	%  errorCovariance == v*d*v'
	%  errorCovariance^0.5 == v*sqrt(d)*v'
	%  errorCovariance^(-0.5) == v*diag(1./sqrt(diag(d)))*v'
	d = diag(d);
	if(any(d<=0))
		error('The covariance matrix is not positive definite');
	end;
	correctionMatrix = v*diag(1./sqrt(d))*v';
	clear v d;
else 
	clear correctionMatrix;
	correctionMatrix = [];
end;

%%	################################# covariates processing #################################

% Add constant as a covariate
cvrt.SetNanRowMean();
cvrt.CombineInOneSlice(); %
cvrt = [ones(1,snps.nCols);cvrt.dataSlices{:}];

% Correct for the error covariance structure
if(~isempty(correctionMatrix))
	status('Rotating cvrt based on the errorCovariance matrix');
	cvrt = cvrt * correctionMatrix;
end;

% Orthonormalize covariates
status('Orthonormalizing covariates');
[q, r] = qr(cvrt',0);
if(min(abs(diag(r))) < eps(class(r))*snps.nCols())
	error('Colinear or zero covariates detected.');
end;
cvrt = q';
clear q;

%%	################################# gene expression processing #################################

% Impute gene expression
status('Imputing missing expression');
gene.SetNanRowMean();
% Correct for the error covariance structure
if(~isempty(correctionMatrix))
	status('Rotating expression based on the errorCovariance matrix');
	gene.RowMatrixMultiply(correctionMatrix);
end;

gene.RowStandardizeCentered();

% Orthogonolize expression w.r.t. the covariates
status('Orthogonolizing expression w.r.t. covariates');
for sl = 1:gene.nSlices
	slice = gene.dataSlices{sl};
	slice = slice - (slice*cvrt')*cvrt;
	gene.dataSlices{sl} = slice;
end;
clear sl d slice;
status('Standardizing expression');
gene.RowRemoveZeroEps();
gene.RowStandardizeCentered();

%%	################################# Prepare for main loop    #################################

nSamples = snps.nCols();
nGenes = gene.nRows;
nSnps  = snps.nRows;
nCov = size(cvrt,1);
% nVarTested = length(snps_list); % set in case(useModel)
% dfNull = nSamples - nCov;
dfFull = nSamples - nCov - nVarTested;

%snameSlices = snps_list.rowNameSlices;
%gnameSlices = gene.rowNameSlices;

switch useModel
	case modelANOVA
		if(pvOutputThreshold >= 1)
			r2Thresh = 0;
		else
			fThresh = finv(1-pvOutputThreshold, nVarTested,dfFull);
			r2Thresh = fThresh * nVarTested ./ (dfFull + fThresh*nVarTested);
			clear fThresh;
		end;
	case modelLINEAR
		if(pvOutputThreshold >= 1)
			rThresh = 0;
		else
			tThresh = -tinv(pvOutputThreshold/2,dfFull);
			rThresh = sqrt(  tThresh.^2 ./  (dfFull + tThresh.^2)  );
			clear tThresh;
		end;
end;

gene_names = vertcat(gene.rowNameSlices{:});
snps_names = vertcat(snps.rowNameSlices{:});
FDR_collection = cell(snps.nSlices,gene.nSlices);
FDR_total_count = 0;

totalCount = nGenes*nSnps;
dumpCount = 0;

%% ################################# Main loop #################################

status('Performing eQTL analysis.');
snps_offset = 0;
tic
for sc = 1:snps.nSlices
	gene_offset = 0;
	
	% prepare snps / dummies
	cursnps = snps_process( snps.dataSlices{sc} );
	for d = 1:length(cursnps)
		if(~isempty(correctionMatrix))
			cursnps{d} = cursnps{d} * correctionMatrix;
		end;
		cursnps{d} = cursnps{d} - (cursnps{d}*cvrt')*cvrt;
		for w = 1:(d-1)
			cursnps{d} = cursnps{d} - ...
				bsxfun(@times, sum(cursnps{d}.*cursnps{w},2), cursnps{w});
		end;
		cursnps{d} = RowStandardizeCentered(cursnps{d});
	end;

	nrcs = size(cursnps{1},1);

	for gc = 1:gene.nSlices
		curgene = gene.dataSlices{gc};
		nrcg = size(curgene,1);
		
		switch useModel
			case modelLINEAR
			cursnps{1}
			curgene
			pause
				cor = (cursnps{1}*curgene');

				cor
				pause
				select = (abs(cor) >= rThresh);
				r = cor(select);
				test = r.*sqrt( dfFull ./ (1-r.^2));
				pv = tcdf(-abs(test),dfFull)*2;
			case modelANOVA
				r2 = (cursnps{1}*curgene').^2;
				for d = 2:nVarTested
					r2 = r2 + (cursnps{d}*curgene').^2;
				end;
				select = (r2 >= r2Thresh);
				rsub = r2(select);
				test = rsub./(1-rsub) * (dfFull/nVarTested);
				pv = fcdf(1./test, dfFull, nVarTested);					
		end;
		[sind gind] = find(select);
		FDR_collection{sc,gc} = [snps_offset+sind, gene_offset+gind, test, pv];
		dumpCount = dumpCount + size(FDR_collection{sc,gc},1); %signifs_sic(1); 

		gene_offset = gene_offset + nrcg;
		FDR_total_count = FDR_total_count + nrcg*nrcs;

		disp([num2str(floor(FDR_total_count/totalCount*1000)/10,'%3.1f') '% done, ' num2str(dumpCount) ' eQTLs found.']);
	end;
	snps_offset = snps_offset + nrcs;
end;
t2 = toc;
clear cor r2 select test pv cursnps curgene

status('Calculating FDR');

FDR_collection = vertcat(FDR_collection{:});

[~, order] = sort(FDR_collection(:,4));
FDR_collection = FDR_collection(order,:);
clear order;
FDR = FDR_collection(:,4)*FDR_total_count./(1:size(FDR_collection,1))';
FDR(end) = min(FDR(end),1);
for i=(size(FDR_collection,1)-1):-1:1
	FDR(i) = min(FDR(i),FDR(i+1));
end;

status('Saving results');
[fid, msg] = fopen(Output_file_name,'w');
if(fid == -1)
	error(msg)
end;
fprintf(fid,'SNP\tgene\t%s\tp-value\tFDR\r\n', statistic_name);
step = 1000;
for part = 1:ceil(size(FDR_collection,1)/step)
	fr = (part-1)*step + 1;
	to = min(part*step,size(FDR_collection,1));

	sdump = snps_names(FDR_collection(fr:to,1));
	gdump = gene_names(FDR_collection(fr:to,2));
	adump = [sdump gdump num2cell([FDR_collection(fr:to,3:4), FDR(fr:to)])]';
	fprintf(fid,'%s\t%s\t%f\t%e\t%e\r\n',adump{:});
end;
fclose(fid);
status('');
clear cor select adump pdump pv testmatrix

disp([num2str(dumpCount) ' eQTLs saved in ' Output_file_name]);
compexity = 2 * nSamples * nGenes * nSnps * nVarTested;
gflops = compexity / t2 / 1e9;

disp(['Performance: ' num2str(gflops) ' gflops']);
end

function status_T(text,~)
	if(nargin==1)
		disp(['Task	finished in ' num2str(toc) ' seconds']);
	end;
	disp(text);
	tic;
end

function status_F(~,~)
end

function y = impute_row_mean(x)
	if(any(x))
		rowmean = nanmean(x,2);
		rowmean(isnan(rowmean)) = 0;
		for j=1:size(x,2)
			where1 = isnan(x(:,j));
			x(where1,j) = rowmean(where1);
		end
	end;
	y = {x};
end

function dummies = snps_split_for_ANOVA(x)
	md = mode(x,2); % Use less frequent values for dummies to reduce colinearity with the constant.
	x(bsxfun(@eq,x,md)) = NaN; 			
	dummies = cell(2,1);
	for d=1:2
		md = mode(x,2);
		dummies{d} = bsxfun(@eq,x,md);
		x(dummies{d}) = NaN;
	end;
	if(~all(isnan(x(:))))
		error('More than 3 SNP values encountered. Not ok for ANOVA as it is coded right now.');
	end;
end

function y = RowStandardizeCentered(x)
	div = sqrt(sum(x.^2,2));
	div(div==0) = 1;
	y = bsxfun(@rdivide, x, div);
end
