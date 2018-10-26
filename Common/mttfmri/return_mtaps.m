function nseq = return_mtaps(baseVal,powerVal);

%  return_mtaps(baseVal,powerVal);
%  
%  returns the number of m-sequences available for a given base and
%  power.
%  
%  Send comments and questions to ttliu@ucsd.edu

if baseVal==2,
  switch powerVal,
   case 2, nseq = 1;
   case 3, nseq = 2;
   case 4, nseq = 2;
   case 5, nseq = 6;
   case 6, nseq = 6;
   case 7, nseq = 18;
   case 8, nseq= 16;
   case 9, nseq = 48;
   case 10, nseq = 60;
	  
   case {11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,39},nseq= 1;
   otherwise error(sprintf('M-sequence %.0f^%.0f is not defined',baseVal,powerVal))
  end;   
elseif baseVal==3,
  switch powerVal,
   case 2, nseq = 2;
   case 3, nseq = 4;
   case 4, nseq = 8;
   case 5, nseq = 21;
   case 6, nseq = 24;
   case 7, nseq = 27;
   otherwise error(sprintf('M-sequence %.0f^%.0f is not defined',baseVal,powerVal))
  end;   
elseif baseVal==5,
  switch powerVal,
   case 2, nseq = 4;
   case 3, nseq = 15;
   case 4, nseq = 48;
   otherwise error(sprintf('M-sequence %.0f^%.0f is not defined',baseVal,powerVal))
  end	
elseif baseVal == 4  | baseVal == 7 | baseVal == 8 | baseVal == 9 | baseVal ...
      == 11
  switch powerVal,
   case 2, if baseVal == 11; nseq = 18; else nseq = 1; end
   case 3, nseq = 1;
   case 4, nseq = 1;
   otherwise error(sprintf('M-sequence %.0f^%.0f is not defined',baseVal,powerVal))
  end	
elseif baseVal == 13
  switch powerVal,
   case 2, nseq = 1;
   case 3, nseq = 1;
   otherwise error(sprintf('M-sequence %.0f^%.0f is not defined',baseVal,powerVal))
  end	  
end;

