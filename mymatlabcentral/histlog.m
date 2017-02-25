function histlog(data1,binnum,plottype,dottype)

if nargin<4, dottype='s'; end
if nargin<3, plottype=1; end
if nargin<2, binnum=20; end
[y,x]=hist(data1,binnum);

%y=y./sum(y);


switch plottype
    case 1
        h=bar(x,y,1.0);
        set(gca,'yscale','log')
        set(h,'BaseValue',0.02)
    case 2
        plot(x,y,dottype);
        set(gca,'xscale','log')
        set(gca,'yscale','log')
    case 3
        
        % http://www.mathworks.com/support/solutions/en/data/1-2ZUTKK/?solution=1-2ZUTKK
        %x = -5:0.1:5;
        %y = randn(100000,1);
        hist(data1,binnum);
        x
        %Workaround
        ph = get(gca,'children');
        vn = get(ph,'Vertices');
        vn(:,2) = vn(:,2) + 1;
        set(ph,'Vertices',vn);
        set(gca,'yscale','log')
    case 4
        hist(data1,binnum);
        x
        %Workaround
        ph = get(gca,'children');
        vn = get(ph,'Vertices');
        vn(:,2) = vn(:,2) + 1;
        set(ph,'Vertices',vn);
        set(gca,'yscale','log')
        set(gca,'xscale','log')
        
end

