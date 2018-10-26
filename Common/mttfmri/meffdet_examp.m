% example script using calc_meffdet.m to 
% look at performance of some  m-sequences and random designs
%
% Send comments and questions to ttliu@ucsd.edu


% specify design parameters
ndesigns = 1;  % number of random designs to try out, 1,000 is a
               % reasonable number, although the more the better
nevents = 10;  % number of trial types, not including control
npts = 121;    % number of reps in the design
nummods = 7;   % number of  points in hemodynamic response function
dt = 3;        % TR of the experimental design in seconds


% make up the hemodynamic response function
t_ax = 0:dt:(nummods-1)*dt;
h = hemoresp(t_ax,1.2,3,1);

% make up matrix to project out nuisance functions
lorder = 3; % degree of Legendre polynomial to use
            % 0 = DC, 1 = DC + Linear, 2 = DC +Lin + Quad, 3 = DC+Lin+Quad+Cubic
S = legendremat(lorder,npts);
PS = eye(npts,npts)-S*pinv(S);

% make up the random designs
randmat = unidrnd(npts,npts,ndesigns);
[ignore,prandmat] = sort(randmat);
basevecmat = rem(prandmat,nevents + 1);

% optional m-sequence generation;
% m-sequences are good for achieving close to
% the theoretical maximum efficiency
mparam.base = 11;mparam.power = 2;mparam.shift = 0;
msmat = gen_mseq(mparam,npts);
basevecmat = [msmat basevecmat];

ndesigns = size(basevecmat,2);

% user defined contrast matrix
% number of columns must be equal to nevents
% number of rows can be anything, i.e. uC need not be full rank

%event  1   2   3   4    5   6   7   8  9   10
uC = [  0   1   0   0    0   0   0   0  0   0;
	0   0   1   0    0   0   0   0  0   0;
	0   0   0   1    0   0   0   0  0   0;
	0   0   0   0    1   0   0   0  0   0;
	0   0   0   0    0   1   0   0  0   0;
	0   0   0   0    0   0   1   0  0   0;	
	0   0   0   0    0   0   0   1  0   0;	
	0   0   0   0    0   0   0   0  1   0;	
	0   0   0   0    0   0   0   0  0   1;	
	0   0   0   0    1/3 1/3 1/3 -1/3 -1/3 -1/3;
	0   0   0   0    1   0   0  -1  0   0;	
	0   0   0   0    0   1   0   0 -1   0;	
	0   0   0   0    0   0   1   0  0  -1;	
	0   1   0   0    0   0   0  -1  0   0;	
	0   1   0   0    0   0   0   0 -1   0;	
	0   0   1   0    0   0   0  -1  0   0;	
	0   0   1   0    0   0   0   0  0  -1;	
	0   0   0   1    0   0   0   0 -1   0;		
	0   0   0   1    0   0   0   0  0  -1;		];
	


% initialize efficiency and power for user-defined contrasts
ueffmat = NaN*ones(ndesigns,1);
udetmat =ueffmat;
% initialize efficiency and power for general contrasts
effmat = NaN*ones(ndesigns,1);
detmat =effmat;



% initialize other variables
V = []; % empty matrix for noise correlation signifies uncorrelated noise
mode = 1; % mode = 1 will calculate both general contrasts and
          % user-defined contrasts
	  
% Theoretical Maximum values, these are valid for
% general contrasts only. M-sequences come very close to the
% theoretical efficiency when lorder = 0.
% Note that the detection power of a random sequence is
% approximately given by: maxdet/nummods
maxeff = npts/(2*(nevents + 1)*nummods); % maximum theoretical efficiency
maxdet = npts*nummods/(2*(nevents + 1)); % maximum theoretical detection power

for k = 1:ndesigns
  k
      [thiseff,thisdet,this_ueff,this_udet] = ...
       calc_meffdet(basevecmat(:,k),nummods,nevents,PS,h,V,uC,mode);
   if mode >= 1
     ueffmat(k) = this_ueff(1);udetmat(k) = this_udet(1);
   end
   if mode <=1
     effmat(k) = thiseff(1);detmat(k) = thisdet(1);
   end
end;


% Now we can print out the best design, based on whatever metric we
% like.
% here we pick the design with the best overall efficiency.
fname = 'besteff.dat';
[m,index] = max(effmat)
printdesign(basevecmat,index,fname);

% here is design with best overall efficiency with user contrasts
fname = 'bestueff.dat';
[m,index] = max(ueffmat)
printdesign(basevecmat,index,fname);

% here is design with best overall detection power with user contrasts
fname = 'bestudet.dat';
[m,index] = max(udetmat)
printdesign(basevecmat,index,fname);

% It's sometimes useful to check a design that is stored in
% a file
disp('checking design')
[e,d,ue,ud] = check_mdesign(fname,nummods,nevents,PS,h,V,uC,mode);

% the two rows should be the same.
[e,d,ue,ud;
 effmat(index) detmat(index) ueffmat(index) udetmat(index)]
