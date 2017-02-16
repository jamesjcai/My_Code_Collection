function wsgi_batch_by_id(varargin)
N=1500;
n=50;
idstart=1:n:N-1;
idend=idstart+n-1;

% http://www.mathworks.com/matlabcentral/answers/53287-how-to-pass-arguments-to-matlab-executables
userArg1 = varargin{1};
% userArg2 = varargin{2};
id=str2double(userArg1);
if isnan(id), error('xxx'); end
if id>length(idstart), error('yyy'); end


for k=idstart(id):idend(id)
    k
end
disp('Done.')



    
