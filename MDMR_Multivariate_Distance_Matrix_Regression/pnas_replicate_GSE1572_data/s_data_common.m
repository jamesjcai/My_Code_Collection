load genelist genelist
load peer_Adipose_Subcutaneous_350_294.mat residuals
load ../smpl_Adipose_Subcutaneous_350
i=vrace==3;
uid=vid(i);
uage=vage(i);
ugender=vgender(i);

d=residuals';
g=genelist;
X=uage;




