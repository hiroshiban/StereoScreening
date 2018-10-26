function ent = calcent2(seq,order);
  
%calculate the conditional entropy of a sequence.
% this function is used to double-check  calcent.m 
% it should only be used when (nevents+1)^order is 
% a manageable number
%  
% appears to be slightly slower than calcent.m

  
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

possible_codes = 0:(b^order-1);
ncodes = length(possible_codes);
pmat = [];
for icode = 1:ncodes;
  code = possible_codes(icode);
  code_inds = find(codes == code);
  code_freq = length(code_inds);
  if code_freq > 0
    for result = 0:nevents
      nresult = sum(results(code_inds) == result);
      if(nresult > 0)
	pmat = [pmat;code result nresult nresult/code_freq];
      end
      
    end
  end
end

jointprob = pmat(:,3)/npos;
condprob = pmat(:,4);
ent = -sum(jointprob.*log2(condprob));

keyboard;






