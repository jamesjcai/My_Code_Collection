function i_halfimage
a=dir('*.png');
[~,idx] = sort([a.datenum]);
a = a(idx);
for k=1:length(a)
    f1=a(k).name;
    A=imread(f1);
    [x,y,~]=size(A);
    if y>3000        
        B=imcrop(A,[0,0,round(y/2),x]);
        imwrite(B,sprintf('img/%d_%s',k,f1));
    end
    k
end
