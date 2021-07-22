pw1=pwd;
if ~exist(fullfile(pw1,'scfft_testdata.mat'),'file')
    try
        disp('Downloading...scfft_testdata.mat');
        websave('scfft_testdata.mat','https://github.com/jamesjcai/jamesjcai.github.io/raw/master/data/scfft_testdata.mat');
    catch
        error('Cannot download scfft_testdata.mat');        
    end
end
load scfft_testdata.mat
X2=sc_impute(X);
%%
% k=200
% [x,idx]=sort(X2(k,:));
for k=200:220
    x=X2(k,:);
    % x=x(idx);
    L=length(x);
    y = fft(x);

    % P2 = abs(y/L);
    P2 = y.*conj(y)/L;
    P1 = P2(1:round(L/2)+1);
    P1(2:end-1) = 2*P1(2:end-1);

    Fs=1000;
    f = Fs*(0:round(L/2))/L;
    
    
    figure;
        subplot(3,1,1)
        plot(x);

    subplot(3,1,2)
        zy=-log10(P1);
        zy(zy<0)=0;
        plot(f,zy)
        yline(5,'r-')
        title('Single-Sided Amplitude Spectrum of X(t)')
        xlabel('f (Hz)')
        ylabel('|P1(f)|')
        xlim([-10 500])
        title(genelist(k))

    subplot(3,1,3)
        title('denoised')
        zy=-log10(2*P2);
        zy(zy<0)=0;
        %idx=zy>5;
        [~,ix]=maxk(zy,2);
        idx=false(size(P2));
        idx(ix)=true;
        if any(idx)
            PSDclean=P2.*idx;
            fhat = idx.*y;
            ffilt = ifft(fhat);
            plot(zscore(x));
            hold on
            plot(zscore(ffilt));
        end    
    
end
