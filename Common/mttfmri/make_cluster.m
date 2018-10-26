% Example of generating many realizations (npaths) of 
% permuted block and clustered m-sequence designs.
% Currently can do either 2 or 4 trial types, see multsim3.m for 
% example of expansion to more trial types
% 
% NOTE: in this example version of the code, npaths = 100 . In the
%       paper, npaths = 1000 (which takes a long time to compute). 
%
% doblock = 1: permuted block designs
% doblock = 0: clustered m-sequences
%
% 

% design parameters
if ~exist('nevents');nevents = 4;end
if ~exist('docorr');docorr = 0;end;
if ~exist('doblock');doblock = 0;end;
npts = 240;
combos = nchoosek(1:nevents,2); % number of combinations of contrasts
ncontrasts = size(combos,1);
neffdet = nevents + ncontrasts + 1;
baseveclen = npts;

if doblock
% make up a simple alternating block design
        numones = npts/(nevents + 1);numblocks = 2;
	sizeofblock = round(numones/numblocks);
	blockspacing = npts/(numblocks);
	basevec = zeros(npts,1);
	for k = 1:nevents
	  thisbasevec =zeros(baseveclen,1);
	  startlocation = 1 + (k-1)*sizeofblock;
	  blocklocs = startlocation:blockspacing:npts;
	  thisbasevec(blocklocs) =1;
	  thisbasevec = conv(thisbasevec,ones(sizeofblock,1));
	  basevec(find(thisbasevec(1:npts))) = k;
	end
	lastnonzero = max(find(basevec > 0));
	nlastzeros = npts-lastnonzero;
	thisoffset = round(nlastzeros/2);
	basevec = shiftvec(basevec,thisoffset);

	ms = basevec;

else % default is m-sequence clumping
% m-sequence generation
if nevents == 4
  mparam.base = 5;mparam.power = 3; % gives length 124 point m-seq
  mparam.shift = 0;
  ms = gen_mseq(mparam,npts);  
elseif nevents == 2
  mparam.base = 3;mparam.power = 5; % gives length 242 m-seq
  mparam.shift = 0;
  ms = gen_mseq(mparam,npts);  
end  

end

% set the random number generator to a known state
load rstate2mult;
rand('state',newrstate);

nummods = 15;
h = hemoresp(0:(nummods-1),1.2,3,1);
if doblock
  nperm_clump = 100;
  npaths = 100;
  nmodes = 1;
else
  nperm_clump = 30;
  npaths = 100;
  nmodes = 1;
end


lorders = [0];
norders = length(lorders);
cleff = NaN*ones(nmodes,norders,nperm_clump,npaths,neffdet); 
cldet = cleff;
clstim =  NaN*ones(nmodes,nperm_clump,npaths,npts);
PS = NaN*ones(norders,npts,npts);
alphavec = [0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8];
alpha = alphavec(1); % just do uncorrelated noise for now. 
V = toeplitz(alpha.^(0:npts-1));
Vi = inv(V);
for iorder = 1:norders
  lorder = lorders(iorder);
  S = legendremat(lorder,npts);
%  PS(iorder,:,:) = eye(npts,npts)-S*pinv(S);
  PS(iorder,:,:) = Vi-Vi*S*inv(S'*Vi*S)*S'*Vi;
end



for mode = 1:nmodes
  for path = 1:npaths
    [mode,path]
    if docorr
      alpha = alphavec(path);
      V = toeplitz(alpha.^(0:npts-1));
      Vi = inv(V);
      %normalize Vi so that trace is npts
      Vi = Vi*npts/trace(Vi);
      S = legendremat(0,npts);
      iorder = 1;
      PS(iorder,:,:) = Vi-Vi*S*inv(S'*Vi*S)*S'*Vi;
    end
    % initialize the starting design
    seq = ms(:,1);
  thistype = 1;
  for k = 1:nperm_clump
    % calculate eff/det for the sequence
    for iorder = 1:norders
     [thiseff,thisdet] = calc_meffdet(seq,nummods,nevents,squeeze(PS(iorder,:,:)),h);
     cleff(mode,iorder,k,path,:) = thiseff;
     cldet(mode,iorder,k,path,:) = thisdet;
    end
    clstim(mode,k,path,:) = seq;
    
    if doblock % de-clump the sequence
      looping = 1;
      while looping % loop until we find two points with different values
	p1 = unidrnd(npts);
	p2 = unidrnd(npts);
	if seq(p1) ~= seq(p2)
	  tmp = seq(p1);
	  seq(p1)=seq(p2);
	  seq(p2)=tmp;
	  looping = 0;
	end
      end
      
    else
    %clump the sequence
      seq = clumpvec(seq,thistype,mode);
    end
    thistype = thistype + 1;
    if thistype > nevents
      thistype = 1;
    end
   end
  end
end


Q = nevents;
npts = 240;
tmaxeff = npts/(2*(Q+1)*nummods);
tmaxdet = npts*nummods/(2*(Q+1));
neff = cleff/tmaxeff;
ndet = cldet/tmaxdet;
if doblock
  pre = 'bl';
  bleff = cleff;
  bldet = cldet;
  blstim = clstim;
  varnames = ['bleff bldet blstim tmaxeff tmaxdet'];
  filename = sprintf('%s%dl%snpe%dnp%dnb%d',pre,nevents, ...
		   no_blanks(num2str(lorders)),nperm_clump,npaths,numblocks);
else
  pre = 'cl';
  varnames = ['cleff cldet clstim tmaxeff tmaxdet'];
  filename = sprintf('%s%dl%snpe%dnp%d',pre,nevents, ...
		   no_blanks(num2str(lorders)),nperm_clump,npaths);
end





eval(sprintf('save %s %s',filename,varnames));
