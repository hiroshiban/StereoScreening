function [ent,ent0] = calcent(seq,order);
  
%calculate the conditional entropy of a sequence.
%
% [ent,ent0] = calcent(seq,order);
%   
%  Inputs
%    seq:  design vector
%  order:  order of conditional entropy
%
%  Outputs
%    ent: conditional entropy
%   ent0: ordinary entropy
%
% send questions to ttliu@ucsd.edu  

seq = seq(:);
nevents = max(seq);  
b = nevents + 1;
seqlen = length(seq);
basevec = b.^((order-1):-1:0);

npos = seqlen-order;

codes = NaN*ones(npos,1);         %conditional codes
results = seq((order+1):seqlen);  %resulting code

for k = 1:npos
  codes(k) = basevec*seq(k+(0:(order-1)));
end

% sort the codes
[scodes,icodes] = sort(codes);
% determine how many distinct codes there are
% in the sequence, this is necessary when nevents and order are large
iresults = results(icodes);% resort the results too
dscodes = diff(scodes);
distinct_codes_ind =  find(dscodes > 0);
span = [1; distinct_codes_ind+1];
distinct_codes = scodes(span);
code_freq = diff([span;(npos+1)]);

bigspan = [span;(npos+1)];

% figure out probabilities of occurrence
pmat = [];
ncodes = length(distinct_codes);
for icode = 1:ncodes
  code = distinct_codes(icode);
  codespan = bigspan(icode):(bigspan(icode+1)-1);
  for result = 0:nevents
    nresult = sum(iresults(codespan) == result);
    if(nresult > 0)
     pmat = [pmat;code result nresult nresult/code_freq(icode)];
    end
  end
end

jointprob = pmat(:,3)/npos;
condprob = pmat(:,4);
ent = -sum(jointprob.*log2(condprob));

% also calculate ordinary entropy
codeprob = NaN*ones(b,1);
for ii = 1:b;
  codeprob(ii) = sum(seq == ii-1);
end;  
codeprob = codeprob/seqlen;
ent0 = -sum(codeprob.*log2(codeprob));





