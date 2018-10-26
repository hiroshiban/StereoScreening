function lf = laguerre(tax,alpha,order);


% generate Laguerre functions with parameter alpha and specified orders
%  lf = laguerre(tax,alpha,order);
%
% send questions to ttliu@ucsd.edu  
  
  norder = length(order);
  lf = NaN*ones(length(tax),norder);
  ind = 1;
  for j = order
    %a1 = alpha.^(tax-j)/2*(1-alpha)^0.5;
    a1 = sqrt(1-alpha^2).*alpha.^(tax-j);
    b1 = 0;
    for k = 0:(j-1)
      b1 = b1 + (-alpha)^k*tchoose(tax+k-1,j-1)*nchoosek(j-1,k)*alpha^k;
    end

    lf(:,ind) = a1(:).*b1(:);
    
    ind = ind + 1;
  end
  


function tc = tchoose(timevec,k);
  tlen = length(timevec);
  tc = NaN*ones(tlen,1);
  for ii = 1:tlen
    if timevec(ii) >= k
     tc(ii) = nchoosek(timevec(ii),k);
    else
     tc(ii) = 0;
    end
  end

  
