function savestrcell2file(X)

fid=fopen('XXX_SAVESTRCELL2_temp.txt','w');
for k=1:length(X)
fprintf(fid,'%s\n',X{k});
end
fclose(fid);
ls *.txt
%dos('notepad XXX_SAVESTRCELL2_temp.txt')
edit XXX_SAVESTRCELL2_temp.txt