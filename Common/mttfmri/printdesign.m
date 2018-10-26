function z = printdesign(basevecmat,index,fname);
  
%    z = printdesign(basevecmat,index,fname);
%
%  utility function to print out a file with the 
%  a selected design.
  
fid = fopen(fname,'w');
if(strcmp(computer,'MAC'))
  fprintf(fid,'%d\r',basevecmat(:,index));
else
  fprintf(fid,'%d\n',basevecmat(:,index));
end
fclose(fid);