function [x,y] = stimpatch(stimulus);


 nstim = length(stimulus);
 x = [0 1 1 0];
 x = repmat(x(:),1,nstim);
 addx = ones(4,1)*(0:(nstim-1));
 x = x + addx;
 y = [ 0 0 1 1];
 y = repmat(y(:),1,nstim);


 
 
 
 
 
 
