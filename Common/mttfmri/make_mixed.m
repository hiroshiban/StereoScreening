% Make up mixed designs
% Currently set up to handle 2, 3 or 4 trial types.
% 

if ~exist('nevents');nevents = 4;end
if ~exist('lorders');lorders = [0 1 2 3];end;

norders = length(lorders);
npts = 240;nummods = 15;
baseveclen = npts;
PS = NaN*ones(norders,baseveclen,baseveclen);
for iorder = 1:norders
  lorder = lorders(iorder);
  S = legendremat(lorder,baseveclen);
  PS(iorder,:,:) = eye(baseveclen,baseveclen)-S*pinv(S);
end


numblocksvec = [ 1 2];
nnblocks = length(numblocksvec);
numblocks = 1;
blklenvec = numblocks*(nevents+1)*(1:fix(npts/(nevents+1)));
nlens = length(blklenvec);
if norders == 1;
neffmix = NaN*ones(nlens,nnblocks);
ndetmix = NaN*ones(nlens,nnblocks);
else
neffmix = NaN*ones(nlens,nnblocks,norders);
ndetmix = NaN*ones(nlens,nnblocks,norders);  
end  
mixedstim = NaN*ones(nlens,npts,nnblocks);
n_entorders = 3;
entmix = NaN*ones(nlens,nnblocks,n_entorders);


switch nevents
 case 2
  mparam.base = 3;mparam.power = 5;
  mparam.shift = 0;
 case 3
  mparam.base = 4;mparam.power = 4;
   mparam.shift = 0;
 case 4
   mparam.base = 5;mparam.power = 3;
   mparam.shift = 0;
end

for iblock = 1:nnblocks
  numblocks = numblocksvec(iblock);
  blklenvec = numblocks*(nevents+1)*(1:fix(npts/(nevents+1)/numblocks));
  nlens = length(blklenvec);
for ilen = 1:nlens
  ilen
  blklen = blklenvec(ilen);
  trunclength = npts-blklen;
  ms = gen_mseq(mparam,trunclength);
  blknpts = blklen;baseveclen = blknpts;
  numones = blklen/(nevents + 1);
  sizeofblock = round(numones/numblocks);
  blockspacing = blknpts/(numblocks);
  basevec = zeros(baseveclen,1);

  for k = 1:nevents
    thisbasevec =zeros(baseveclen,1);
    startlocation = 1 + (k-1)*sizeofblock;
    blocklocs = startlocation:blockspacing:blknpts;
    thisbasevec(blocklocs) =1;
    thisbasevec = conv(thisbasevec,ones(sizeofblock,1));
    basevec(find(thisbasevec)) = k;
  end
%  basevec = circshift(basevec,-10);

seq = [ms(:,1);basevec];
for iorder = 1:norders;
  thisPS = squeeze(PS(iorder,:,:));
  [eff,det] = calc_meffdet(seq,nummods,nevents,thisPS);
  
if nevents == 1;
  neffmix(ilen,iblock) = eff(1);
  ndetmix(ilen,iblock) = det(1);
else
  
  neffmix(ilen,iblock,iorder) = eff(1);
  ndetmix(ilen,iblock,iorder) = det(1);
end


  %keyboard; pause;
end
mixedstim(ilen,:,iblock) = seq(:);
end

%calculate entropies of designs
for entorder = 1:n_entorders
  entmix(1:nlens,iblock,entorder) = calcentvec(squeeze(mixedstim(1:nlens,:,iblock)),entorder);
end
end % numblocks loop

Q = nevents;
tmaxeff = npts/(2*(Q+1)*nummods);
tmaxdet = npts*nummods/(2*(Q+1));
neffmix = neffmix/tmaxeff;
ndetmix = ndetmix/tmaxdet; 



if norders == 1
  filename = sprintf('mixed_ne%dnp%d',nevents,npts);
  eval(sprintf('save %s neffmix ndetmix nummods nevents npts mixedstim entmix',filename));
else
  filename = sprintf('mixed_ne%dnp%dno%d',nevents,npts,norders);
  eval(sprintf('save %s neffmix ndetmix nummods nevents npts mixedstim entmix',filename));  
end  

