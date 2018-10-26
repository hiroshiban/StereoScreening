% batch file to run simulations for looking at 
% the effect of basis function expansions. 

doclump = 1; % do clustering of m-sequences
nevents = 2; % Q = 2
numperm = 200;      % number of block permutations
numpermrand = 200;  % number of random sequences to try
nperm_clump = 200;  % number of m-sequence permutations.
basistype = 'friston';
dobasis = 0;
multsim3;
dobasis = 1;
multsim3;


