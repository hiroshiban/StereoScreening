function entvec = calcentvec(seqmat,order);

%  entvec = calcentvec(seqmat,order);
%  
%  calculate conditional entropies of a matrix of sequences
%
%  Inputs:
%    seqmat: matrix of squences, with each row corresponding to a
%             sequence
%    order:  conditional entropy order
%
%  Outputs:
%  
%    entvec: vector of conditional entropies  
%
%  Send questions to ttliu@ucsd.edu
  
  
if size(seqmat,2)> 1
nseq = size(seqmat,1);
else nseq = 1
end  

entvec = NaN*ones(nseq,1);

if nseq > 1
  for iseq = 1:nseq
   entvec(iseq) = calcent(seqmat(iseq,:),order);
  end
else
  entvec = calcent(seqmat,order);
end