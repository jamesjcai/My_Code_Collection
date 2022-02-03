% define figure
f1=figure
ax1=axes('position',[.05 .6 .2 .325]);
ax2=axes('position',[.3 .6 .2 .325]);
ax3=axes('position',[.05 .075 .45 .35]);
ax4=axes('position',[.525 -.20 .475 1.4]);
set(gcf,'Position',[100 100 1400 650])

% define movie
if(movieFLAG)    
    writerObj = VideoWriter([filename,'_Torus.mp4'],'MPEG-4');
    open(writerObj);
end

% define torus grid
r1 = 2;
r2 = 1;
[T1,T2] = meshgrid(0:2*pi/n:2*pi,0:2*pi/n:2*pi);
R = r1 + r2*cos(T2);
Z = r2*sin(T2);
X = cos(T1).*R;
Y = sin(T1).*R;

% define FFT domain grid
x1 = 1:n;
y1 = x1;
x2 = 1:16;
y2 = x2;

for i=1:size(XDAT,2)
    xtilde = reshape(XDATtildeNoise(:,i),n,n);
    x = reshape(XDAT(:,i),n,n);
    set(gcf,'CurrentAxes',ax1)
    imagesc(real(xtilde'),[-1 1]), box on;
    axis xy
    set(gca,'XTick',[1 32 64 96 128],'YTick',[1 32 64 96 128]);
    rectangle('Position',[1 1 16 16],'LineWidth',2)
    set(gcf,'CurrentAxes',ax2)
    imagesc(x2,y2,real(xtilde(x2,y2)'),[-1 1]), box on;  % added transpose!
    axis xy
    set(gca,'XTick',[1 4 8 12 16],'YTick',[1 4 8 12 16]);
    set(gcf,'CurrentAxes',ax3)
    imagesc(abs(x)), box on;
    axis xy
    set(gca,'XTick',[1 32 64 96 128],'YTick',[1 32 64 96 128]);
    set(gcf,'CurrentAxes',ax4)   
    surf(X,Y,Z,abs(x),'EdgeColor','none')
    axis equal
    set(gcf,'Color',[1 1 1]);
    view(100,32)
    colormap(jet);
    drawnow    
    if(movieFLAG)
        frame = getframe(f1);
        writeVideo(writerObj,frame);
    end
end
if(movieFLAG) close(writerObj);
end