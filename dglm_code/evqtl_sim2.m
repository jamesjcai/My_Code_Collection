

N=10000;
testn=200:50:1000;
%testn=300:50:600
naa=1600; nAa=4800; nAA=3600;

T=length(testn);
pv1=nan(T,100);
pv2=nan(T,100);

for t=1:T
    k=1; 
    while k<=100
        y=[1*randn(naa,1); 2+sqrt(2)*randn(nAa,1); 4+sqrt(3)*randn(nAA,1)]+100;
        x=[zeros(naa,1); ones(nAa,1); 2*ones(nAA,1)];
        
        n=testn(t);
        idx=randperm(N);
        i=idx(1:n);
        y=y(i);
        x=x(i);

        %%
        dglmok=false;
        saveR('Data.r','x','y')
        [cstat,p1] = system([getRscript, ' run_dglm.r']);
        if cstat==0
            p1=str2double(p1);
            if ~isnan(p1)
                dglmok=true;
                pv1(t,k)=p1;
            end
        end
        
        %%
        %boxplot(y1,x) % evqtlplot(y1,x,1)
        %{
        mdl = LinearModel.fit(x,y1);
        y2=mdl.Residuals.Raw;
        mdl2 = LinearModel.fit(x,y2.^2);
        mdl2.Coefficients.Estimate(2)    % snp_eff_dispertion
        mdl2.Coefficients.pValue(2)      % snp_dispertion p-value =0.8328
        %}

        %%
        if dglmok
            [b]=regress(y,[ones(size(x)) x]);
            y2=y-(b(1)+b(2)*x);
            [effs, ~, ~, ~, Stats]=regress(y2.^2,[ones(size(x)) x]);
            p2=Stats(3);
            %effsize=effs(2)
            pv2(t,k)=p2;            
            k=k+1;
        end
        fprintf('%d...%d\n',testn(t),k);
    end    
end
readme='sample size = 200:50:1000; pv1 = dglm; pv2 = variABEL';
% save evqtl_sim_res pv1 pv2 readme testn


%%


load evqtl_sim_res pv1 pv2 readme testn
C=sum(pv1<0.05/100,2)./100;
D=sum(pv2<0.05/100,2)./100;


h=figure;
plot(testn,C','o');
hold on
p=polyfit(testn,C',3);
xt=min(testn):max(testn);
f=polyval(p,xt);
%f(f>0.999)=1;
plot(xt,f,'b-');

plot(testn,D','ro');
p=polyfit(testn,D',3);
xt=min(testn):max(testn);
f=polyval(p,xt);
%f(f>0.999)=1;
plot(xt,f,'r-');

pubgraph(h)
yy=ylim;
ylim([yy(1) 1.03])
xlabel('Sample size')
ylabel('Power')

%%
T=3
for t=3:2:10
    %subplot(4,2,t)
    h=figure;
    scatter(-log10(pv2(t,:)),-log10(pv1(t,:)))
    xlabel('variABEL')
    ylabel('dglm')
    vline(-log10(0.05/100))
    hline(-log10(0.05/100))
    refline(1,0);
    title(sprintf('size=%d',testn(t)))
    %pubgraph(h)
    box on
end


