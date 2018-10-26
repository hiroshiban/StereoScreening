function p = optprob(Q,k);
  
%  p = optprob(Q,k);
%
%  Returns optimum frequency of occurrence (p) 
%  for an experiment with Q trial types
%  and contrast weighting parameter 0 <= k <= 1.
%  If only contrasts are of interest: k = 0
%  If only individual triasl are of interest: k = 1  
%
%  examples:
%   
%     optprob(2,0) = 0.5
%     optprob(2,1) = .2929  
  
  p = ( Q*(1-2*k)+Q^2*(k-1)+sqrt(k)*sqrt(Q*(2*k-1)+Q^2*(1-k)))/...
      (Q*(Q-1)*(k*Q-Q-k));