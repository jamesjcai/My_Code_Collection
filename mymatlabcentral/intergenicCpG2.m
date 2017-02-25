function [x,CPG,GC]=intergenicCpG2

CPG=[];
GC=[];
chromset=[1,10,11,12,13,14,15,16,17,18,19,2,20,21,22,3,4,5,6,7,8,9,23,24];


for (chridx=1:24),
    chrid=chromset(chridx);
    fprintf('Chromosome %d\n',chrid);
    T=regioncoord(chrid,'alu',0);
    regions=min(T(:,1)):max(T(:,1));
    n=length(regions);
       
       clear chrmm;
       mmFilename=['C:/biodata/hgenome/Homo_sapiens.NCBI36.40.dna.chromosome.',int2str(chrid),'.mm'];
       chrmm = memmapfile(mmFilename, 'format', 'uint8');
    
    for (k=1:n),
    region=regions(k);
    subT=T(find(T(:,1)==region),:);
    
    totalbp=zeros(1,4);
    totalcpg=zeros(1,6);
    for(j=1:size(subT,1)),
      startn=max(1,subT(j,3)); endn=max(1,subT(j,4));
	  seq=chrmm.Data(startn:endn);
      totalbp=totalbp+[sum(seq==1),sum(seq==2),sum(seq==3),sum(seq==4)];
      numx = i_cpgfreq(seq');
	  totalcpg = totalcpg + numx;
    end    
    CPG=[CPG;totalcpg];
	if (sum(totalbp)>0)
	    GC=[GC;sum(totalbp([2 3]))/sum(totalbp)];
    else
        GC=[GC;0];
    end
    
    end 

end
x=[GC,CPG];




function [num] = i_cpgfreq(s)
	s(find(s>4|s<1))=[];    %remove gaps
	%[Kf] = karlinsig(s);
	%rho=Kf(2,3);
	a=dinuclise16(s);
	b=dinuclise16(s(2:end));
	c=[a,b];

dinucset=[2 3;3 2;4 3;1 3;3 1;3 4];   % CG num TG AC CA GT
m=length(dinucset);
num=zeros(1,m);

for (k=1:m),
      num(k)=sum(c==dinuclise16(dinucset(k,:)));
end
