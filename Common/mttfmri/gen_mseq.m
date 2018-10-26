function msmat = gen_mseq(mparam,mlen,taps);
  
  
% Usage:
%
%
%   msmat = gen_mseq(mparam,mlen,taps);
%  
% Inputs:
%
%  mparam: data structure with
%           mparam.base, mparam.power, mparam.shift
%  mlen  : optional input for desired length of the m-sequence, default = base^power-1
%
%  taps  : optional input for specifying what taps to use, useful
%          for searching for m-sequence tap combos
%
% Outputs:
%  
%  msmat:   matrix with m-sequences, one per column
%
%  Copyright (c) 2002 UC Regents
%
% Send comments and questions to ttliu@ucsd.edu
  
  
if ~exist('taps')
  taps = [];
end
nmseq = return_mtaps(mparam.base,mparam.power);

baselen = mparam.base^mparam.power -1;
if ~exist('mlen');
  mlen = baselen;
end
ms = NaN*ones(mlen,nmseq);
mrepeat = ceil(mlen/baselen);

switch(mparam.base)
 case{2,3,4,5,7,8,9,11,13}
   % do nothing
 otherwise
  error(sprintf('mparam.base ^%d not supported',mparam.base));
end

for ws = 1:nmseq
  ms=mseq2(mparam.base,mparam.power,mparam.shift,ws,taps)';
  ms = ms(:);
  ms = repmat(ms,mrepeat,1);
  lms = length(ms);
  %zeropad or truncate sequence
  if(mlen > lms) % zeropad
    ms = [ms;zeros(mlen-lms,1)];
  else %truncate
    ms = ms(1:mlen);
  end
  msmat(:,ws) = ms;
end