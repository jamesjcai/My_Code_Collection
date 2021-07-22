% https://www.bogotobogo.com/Matlab/Matlab_Tutorial_DFT_Discrete_Fourier_Transform.php
% https://www.youtube.com/watch?v=c249W6uc7ho

x = [2 3 -1 4];
N = length(x);
X = zeros(4,1)
for k = 0:N-1
    for n = 0:N-1
        X(k+1) = X(k+1) + x(n+1)*exp(-j*pi/2*n*k)
    end
end

t = 0:N-1
subplot(311)
stem(t,x);
xlabel('Time (s)');
ylabel('Amplitude');
title('Time domain - Input sequence')

subplot(312)
stem(t,X)
xlabel('Frequency');
ylabel('|X(k)|');
title('Frequency domain - Magnitude response')

subplot(313)
stem(t,angle(X))
xlabel('Frequency');
ylabel('Phase');
title('Frequency domain - Phase response')

X         % to check |X(k)|
angle(X)  % to check phase


%%        


N=4;
x=[2 3 -1 4];

t=0:N-1;
subplot(311)
stem(t,x);
xlabel('Time (s)');
ylabel('Amplitude');
title('Input sequence')

subplot(312); 
stem(0:N-1,abs(fft(x)));  
xlabel('Frequency');
ylabel('|X(k)|');
title('Magnitude Response'); 

subplot(313); 
stem(0:N-1,angle(fft(x)));
xlabel('Frequency');
ylabel('Phase');
title('Phase Response'); 


% https://www.google.com/url?sa=t&rct=j&q=&esrc=s&source=web&cd=&ved=2ahUKEwjUqcD6zODxAhU4EVkFHZu7A5UQFnoECAgQAA&url=https%3A%2F%2Fpeer.asee.org%2Fmatlab-exercises-to-explain-discrete-fourier-transforms.pdf&usg=AOvVaw3NFjx6MlN6_rHC4y7x93v8
% https://www.researchgate.net/post/Can-anyone-provide-me-with-the-MATlab-code-for-DFT-not-the-built-in-function-fft
[Xk] = dft(x)
[Y] = fDFT(x)
[Z]=DFT_nik(x)


function [Xk] = dft(xn)
    xn=xn(:);
    N=length(xn);
    n = 0:1:N-1; % row vector for n
    k = 0:1:N-1; % row vecor for k
    WN = exp(-1j*2*pi/N); % Twiddle factor (w)
    nk = n'*k; % creates a N by N matrix of nk values
    WNnk = WN .^ nk; % DFT matrix
    size(xn)
    Xk = (WNnk*xn);
end


function [Y] = fDFT(x)
%The discrete Fourier transform of input data x
Lx = length(x);
NDFT = 2^nextpow2(Lx); % Next power of 2 from length of x
X =[x zeros(1,NDFT-Lx)];
Y = zeros(1,NDFT);
for k = 0:NDFT-1
Y(k+1) = 0;
for n = 0:NDFT-1
Y(k+1)=Y(k+1)+(X(n+1)*exp((-1j)*2*pi*k*n/Lx));
end
end
end

% http://fweb.wallawalla.edu/class-wiki/index.php/DFT_example_using_MATLAB_-_HW11


function [X]=DFT_nik(x)
%x= input sequence
%X= dft of sequence
% Example x=[1 2 4  1 2 1 2 2 3 2 3];
%If you have any problem or feedback please contact me @
%%===============================================
% NIKESH BAJAJ
% Asst. Prof., Lovely Professional University, India
% Almameter: Aligarh Muslim University, India
% +919915522564, bajaj.nikkey@gmail.com
%%===============================================
N=length(x);
for k=1:N
    X(k)=0;
    for n=1:N
        X(k)=X(k)+x(n).*exp(-1j.*2.*pi.*(n-1).*(k-1)./N);
    end
end

f=0:N-1;
subplot(311)
stem(0:N-1,x)
title('Sequence x (in time domain)')
xlabel('time')
ylabel('Amplitude')
grid;

subplot(323)
stem(f,abs(X))
title('Magnitude of Fourier Coeffients using function')
ylabel('|X|')
grid;

subplot(325)
stem(f,angle(X))
title('Angle of Fourier Coeffients using function')
xlabel('Frquency coefficients')
ylabel('<X')
grid;

subplot(324)
stem(f,abs(fft(x)))
title('Magnitude of Fourier Coeffients using fft function')
ylabel('|X|')
grid;

subplot(326)
stem(f,angle(fft(x)))
title('Angle of Fourier Coeffients using fft function')
xlabel('Frquency coefficients')
ylabel('<X')
grid;
end

% https://github.com/Nikeshbajaj/Discrete-Fourier-Transform/blob/master/DFT_nik.m

% https://www.algorithm-archive.org/contents/cooley_tukey/cooley_tukey.html