for k=1:length(A_pastespecial)
a=A_pastespecial{k};
[b,c,d]=genenamesearchX(a);
fprintf('%s\t%s\t%s\t%s\n',a,b,c,d);
end