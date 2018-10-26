function y = no_blanks(str);

   x = abs(str);
   not_blank = find(x ~= 32);
   y = str(not_blank); 
