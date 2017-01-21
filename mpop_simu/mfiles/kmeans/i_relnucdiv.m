function [d]=i_relnucdiv(seq1,seq2)

d1=i_quickpi(seq1,1);
d2=i_quickpi(seq2,1);
d3=i_quickpi2(seq1,seq2,1);
d=d3./max([d1 d2]);


function [dv,d] = i_quickpi2(seq1,seq2,jc)
    [n1,m1]=size(seq1);
    [n2,m2]=size(seq2);
    
	dv=0; d=0;
	if (n1==1||n2==1), return; end
	  ign=0;  % ignored cell where JC distance cannot be computed
	for i=1:n1
	for j=1:n2
	    p=sum(seq1(i,:)~=seq2(j,:))./m1;
	    if (jc), 
	      if p>=0.75;
            p=0;
            ign=ign+1;
	      else
      		p=(-3/4)*log(1-4*p./3);
	      end
	    end
	    d=d+p;
	end
    end
    n=n1*n2;
	dv=d/(n*(n-1)/2-ign);


function [dv,d] = i_quickpi(seq,jc)
    [n,m]=size(seq);
	dv=0; d=0;
	if (n==1), return; end
	  ign=0;  % ignored cell where JC distance cannot be computed
	for i=1:n-1
	for j=i+1:n
	    p=sum(seq(i,:)~=seq(j,:))./m;
	    if (jc), 
	      if p>=0.75;
            p=0;
            ign=ign+1;
	      else
      		p=(-3/4)*log(1-4*p./3);
	      end
	    end
	    d=d+p;
	end
	end		
	dv=d/(n*(n-1)/2-ign);
    