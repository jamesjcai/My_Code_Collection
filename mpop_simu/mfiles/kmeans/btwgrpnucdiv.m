function [d]=btwgrpnucdiv(hap)

[hap,idx]=i_kmeanssort(hap,2);
[d]=i_relnucdiv(hap(1:idx,:),hap(idx+1:end,:));