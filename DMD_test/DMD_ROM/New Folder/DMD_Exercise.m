close all
clear
clc

%------------- Read the image data -------------------

NumberImages = 100;
NumberDMDModes = 30;
X = [];
I = [];

%------------------ Display data ----------------------

figure;
h = imshow( I );
axis tight;
title('Loading image data...');
disp('Loading image data...');

%------------------- Read data ------------------------

for i=1:NumberImages
    FileName = sprintf('Test%03d.png',i);
    
    I = imread( FileName );
    I = double( I ) / 255;
    
    X = [ X , I(:) ];
    
    set( h , 'CData' , I );
    drawnow;    
    pause(0.01);
end

title('Calculating DMD...');
disp('Calculating DMD ...');
pause(0.01);

%------------------ Data formatting ------------------

ISize = size( I );

X1 = X( : , 1:(end-1) );
X2 = X( : , 2:end );

%----------------------- DMD -------------------------

[ EigenVector , EigenValue ] = DMD( X1 , X2 , NumberDMDModes );

%------------------ Plot DMD modes -------------------

f = figure;
f.WindowState = 'maximized';
k = 1;

for i=1:3
    for j=1:10
        subplot(3,10,k);
        I = reshape( EigenVector(:,k) , ISize );
        imagesc( abs(I) );
        
        D = EigenValue(k,k);
        title( sprintf('EV = %0.2f + i %0.2f', real(D) , imag(D)) );
        
        k = k+1;
    end
end

sgtitle('First 30 DMD with eigen values');
