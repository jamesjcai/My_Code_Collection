classdef SlicedData < handle
	properties
		dataSlices = {};
		rowNameSlices = {};
		columnNames = {};

		fileDelimiter = char(9);
		fileSkipColumns = 1;
		fileSkipRows = 1;
		fileSliceSize = 10000;
		fileOmitCharacters = 'NAna?';
	end % properties
	properties (Access=private)
		fileOmitCharactersReplaceBy = ' ';
	end;
	properties (SetAccess = private, Dependent)
		nRows;
		nCols;
		nSlices;
	end % properties (SetAccess = private)
	methods
		function set.fileDelimiter(obj, delimiter)
			obj.fileDelimiter = delimiter;
			if(obj.fileDelimiter == ' ')
				obj.fileOmitCharactersReplaceBy = char(9);
			else
				obj.fileOmitCharactersReplaceBy = ' ';
			end;
		end;
		function nRows = get.nRows(obj)
			nRows  = sum( cellfun(@(x)size(x,1),obj.dataSlices)); %snps.nRows;
		end % get.nRows
		function nCols = get.nCols(obj)
			if(obj.nSlices == 0)
				nCols = 0;
			else
				nCols = size(obj.dataSlices{1},2); %snps.nRows;
			end;
		end % get.nRows
		function nSlices = get.nSlices(obj)
			nSlices  = length(obj.dataSlices); %snps.nRows;
		end % get.nRows
		function Clear(obj)
			obj.dataSlices = {};
			obj.rowNameSlices = {};
			obj.columnNames = {};
		end
		function RowCenter(obj)
			for sl=1:obj.nSlices
				obj.dataSlices{sl} = bsxfun(@minus, obj.dataSlices{sl}, mean(obj.dataSlices{sl},2));
			end;
		end;
		function RowCenterNan(obj)
			for sl=1:obj.nSlices
				obj.dataSlices{sl} = bsxfun(@minus, obj.dataSlices{sl}, nanmean(obj.dataSlices{sl},2));
			end;
		end;
		function SetNanZero(obj)
			for sl=1:obj.nSlices
				obj.dataSlices{sl}(isnan(obj.dataSlices{sl})) = 0;
			end;
		end;
		function SetNanRowMean(obj)
			for sl=1:obj.nSlices
				select = isnan(obj.dataSlices{sl});
				if(any(select(:)))
					rowmean = nanmean(obj.dataSlices{sl},2);
					rowmean(isnan(rowmean)) = 0;
					for j=1:size(obj.dataSlices{sl},2)
						where1 = isnan(obj.dataSlices{sl}(:,j));
						obj.dataSlices{sl}(where1,j) = rowmean(where1);
					end
				end;
			end;
		end;
		function RowStandardizeCentered(obj)
			for sl=1:obj.nSlices
				div = sqrt(sum(obj.dataSlices{sl}.^2,2));
				div(div==0) = 1;
				obj.dataSlices{sl} = bsxfun(@rdivide, obj.dataSlices{sl}, div);
			end;
		end;
		function RowMatrixMultiply(obj, multiplier)
			for sl=1:obj.nSlices
				obj.dataSlices{sl} = obj.dataSlices{sl}*multiplier;
			end;
		end;
		function RowStandardize(obj)
			obj.RowCenter();
			obj.RowStandardizeCentered();
		end;
		function ColumnSubsample(obj,selection)
			for sl=1:obj.nSlices
				obj.dataSlices{sl} = obj.dataSlices{sl}(:,selection);
			end;
			obj.columnNames = obj.columnNames(selection);
		end;
		function RowSubsample(obj,selection)
			if(obj.nSlices>1)
				error('Sliced matrices are not supported yet. Use CombineInOneSlice first.');
			end;
			obj.dataSlices{1} = obj.dataSlices{1}(selection,:);
			obj.rowNameSlices{1} = obj.rowNameSlices{1}(selection);
		end;
		function CombineInOneSlice(obj)
			if(obj.nSlices>1)
				obj.dataSlices = {vertcat(obj.dataSlices{:})};
				obj.rowNameSlices = {vertcat(obj.rowNameSlices{:})};
			end;
		end;
		function Transpose(obj)
			if(obj.nSlices>1)
				error('Transpose of sliced matrix is not supported yet. Use CombineInOneSlice first.');
			end;			
			obj.dataSlices{1} = obj.dataSlices{1}';
			temp = obj.rowNameSlices{1};
			obj.rowNameSlices{1} = obj.columnNames;
			obj.columnNames = temp;
		end;
		function ResliceCombined(obj)
			if(obj.nSlices>1)
				error('Reslice of sliced matrix is not supported yet. Use CombineInOneSlice first.');
			end;			
% 			oldData = obj.dataSlices{1};
% 			oldRowNames = obj.rowNameSlices{1};
% 			nRows = obj.nRows;
% 			newNSlices = floor((nRows + obj.fileSliceSize -1)/obj.fileSliceSize);
% 			obj.dataSlices = cell(newNSlices,1);
% 			obj.rowNameSlices = cell(newNSlices,1);
% 			for sl=1:newNSlices
% 				obj.dataSlices{sl} = oldData(1 + (sl-1)*obj.fileSliceSize : min(sl*obj.fileSliceSize,end),:);
% 				obj.rowNameSlices{sl} = oldRowNames(1 + (sl-1)*obj.fileSliceSize : min(sl*obj.fileSliceSize,end));
% 			end;
			nRows = obj.nRows;
			newNSlices = floor((nRows + obj.fileSliceSize -1)/obj.fileSliceSize);
			obj.dataSlices    = mat2cell(obj.dataSlices{1}   ,  [obj.fileSliceSize*ones(newNSlices-1,1); nRows - obj.fileSliceSize*(newNSlices-1)],size(obj.dataSlices{1}   ,2));
			obj.rowNameSlices = mat2cell(obj.rowNameSlices{1},  [obj.fileSliceSize*ones(newNSlices-1,1); nRows - obj.fileSliceSize*(newNSlices-1)],size(obj.rowNameSlices{1},2));
		end;
		function copy = Clone(obj)
			copy = SlicedData;

			copy.dataSlices = obj.dataSlices;
			copy.rowNameSlices = obj.rowNameSlices;
			copy.columnNames = obj.columnNames;

			copy.fileDelimiter = obj.fileDelimiter;
			copy.fileSkipColumns = obj.fileSkipColumns;
			copy.fileSkipRows = obj.fileSkipRows;
			copy.fileSliceSize = obj.fileSliceSize;
			copy.fileOmitCharacters = obj.fileOmitCharacters;
		end;
		function SaveText(obj, filename)%, numformat)
			[fid, msg] = fopen(filename,'w');
			if(fid == -1)
				disp(['Error saving to file: ' filename]);
				disp(msg);
				return;
% 			else
% 				msg = '';
			end;
			disp(['Saving file: ' filename]);
% 			if(nargin<=2)
% 				numformat = '%d';
% 			end;
			fprintf(fid,'id');
			fprintf(fid,'\t%s',obj.columnNames{:});
			fprintf(fid,'\n');
			SavedRows = 0;
			for s = 1:obj.nSlices
				for i = 1:size(obj.dataSlices{s},1)
					fprintf(fid,'%s',obj.rowNameSlices{s}{i});
					fprintf(fid,'\t%d',obj.dataSlices{s}(i,:));
					fprintf(fid,'\n');
				end;
				SavedRows = SavedRows + size(obj.dataSlices{s},1);
				disp(['Saved rows: ' num2str(SavedRows) ' of ' num2str(obj.nRows)]);
			end;
			fclose(fid);
			disp('Done saving');
		end;
		function LoadFile(obj, filename)
			%%
			[fid,msg] = fopen(filename,'r');
			if(fid==-1)
				disp(['Error opening file: ' filename]);
				disp(msg);
				return;
			end;
			clear msg;

			for i=1:obj.fileSkipRows
				line = fgetl(fid);
			end;

			if(obj.fileSkipRows == 0)
				line = fgetl(fid);
				fseek(fid,0,'bof');
			end

			positions = find(line == obj.fileDelimiter);
			positions2 = [0,positions,length(line)+1];
			nColumns = length(positions) + 1 - obj.fileSkipColumns;
			obj.columnNames = cell(nColumns,1);
			for i=1:nColumns
				obj.columnNames{i} = line(positions2(obj.fileSkipColumns+i)+1:positions2(obj.fileSkipColumns+i+1)-1);
			end;
			clear positions positions2 i

			%%

			% initial length is 1
			obj.dataSlices = cell(1,1);
			obj.rowNameSlices = cell(1,1);

			curSliceId = 0; % will be incremented soon

			while(true)
				% stretch dataSlices and rowNameSlices if needed
				if(length(obj.dataSlices) <= curSliceId)
					obj.dataSlices{2*curSliceId,1} = [];
					obj.rowNameSlices{2*curSliceId,1} = [];
				end;
				curSliceId = curSliceId + 1;

				% values and names of the current slice
				curRowNames = cell(obj.fileSliceSize,1);
				curRowVals = cell(obj.fileSliceSize,1);

				for i=1:obj.fileSliceSize
					line = fgetl(fid);
					while(ischar(line) && isempty(strtrim(line)))
						line = fgetl(fid);
					end;

					if(line==-1) % end of file
						if(i == 1) % message that the slice is empty
							curRowNames = [];
							curRowVals = cell(0);
						else	   % shrink the slice
							curRowNames = curRowNames(1:(i-1));
							curRowVals = curRowVals(1:(i-1));
						end;
						break;
					end;

					if(nColumns == 0)
						curRowNames{i} = line;
						curRowVals{i} = zeros(1,0);
						continue;
					end;
					
					if(obj.fileSkipColumns > 0)
						position = find(line == obj.fileDelimiter, obj.fileSkipColumns);
						curRowNames{i} = line(1:position(1)-1);
						%position = position(end);
						remainder = line(position(end)+1:end);
					else
						remainder = line;
					end;

					% replace NA with ' '
					for lv = 1:length(obj.fileOmitCharacters)
						remainder(remainder==obj.fileOmitCharacters(lv)) = obj.fileOmitCharactersReplaceBy;
					end;

					nums = textscan(remainder,'%f','Delimiter',obj.fileDelimiter);
					if(length(nums{1}) == nColumns)
						curRowVals{i} = nums{1}';
					else
						curRowVals{i} = [nums{1}', NaN(1,nColumns-length(nums{1}))];
						if(nColumns-length(nums{1}) > 2)
							disp('Parsed line of input file is too short.');
							disp('Problem in line:');
							disp(remainder);
							error('Check "fileOmitCharacters", "fileSkipRows", and "fileSkipColumns"');
						end;
					end;
				end;


				% save the slice in the array
				obj.rowNameSlices{curSliceId} = curRowNames;
				candidate = cell2mat(curRowVals);
% 				candidate8 = int8(candidate);
% 				if(all(candidate(:) == candidate8(:)))
% 					candidate = candidate8;
% 				end;
				obj.dataSlices{curSliceId} = candidate;

				if(isempty(curRowNames))	% there were no slice to read
					curSliceId = curSliceId - 1;
				end;
				% is it is not full, return.
				if(length(curRowNames)<obj.fileSliceSize)
					disp(['Finished reading: ' num2str(obj.nRows()) ]);
% 					disp(['Finished reading: ' num2str((curSliceId-1)*obj.fileSliceSize+length(obj.rowNameSlices{curSliceId})) ]);
					break;
				end;
				disp(['Rows read: ' num2str(curSliceId*obj.fileSliceSize) ]);
			end;
			fclose(fid);
			clear fid;

			obj.rowNameSlices = obj.rowNameSlices(1:curSliceId);
			obj.dataSlices = obj.dataSlices(1:curSliceId);
			
			if(obj.fileSkipRows == 0)
				obj.columnNames = cellstr(num2str(((1:nColumns)'),['Col_%0' num2str(ceil(log10(nColumns+1))) 'd']));
			end;
			if(obj.fileSkipColumns == 0)
				cnt = 0;
				for sl = 1:length(obj.dataSlices)
					nr = size(obj.dataSlices{sl},1);
					obj.rowNameSlices{sl} = cellstr(num2str((cnt+(1:nr)'),['Row_%0' num2str(ceil(log10(obj.nRows()+1))) 'd']));
					cnt = cnt + nr;
				end;
			end;
		end;
		function RowRemoveZeroEps(obj)
			for sl=1:obj.nSlices
				slice = obj.dataSlices{sl};
				amean = mean(abs(slice),2);
				remove = (amean < eps(class(slice))*obj.nCols());
				if(any(remove))
					obj.rowNameSlices{sl} = obj.rowNameSlices{sl}(~remove);
					obj.dataSlices{sl} = slice(~remove, :);
				end;
			end;
		end;
	end% methods
end% classdef