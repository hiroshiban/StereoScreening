function [ne,nd,permind] = compress_clump(nec,ndc);
  
% 030513 - modifying this to spit out permutation indices  
nec = nec(:);ndc = ndc(:);
nbins = 150;
[n,x] = hist(nec,nbins);
dx2 = (x(2)-x(1))/2;
nd = NaN*ones(nbins,3);
ne = nd;
permind = nd;
for k = 1:nbins;
  points = find((nec >= x(k)-dx2) & (nec < x(k)+dx2)); 
  [dmax,ii] = max(ndc(points));

  if ~isempty(dmax)
    nd(k,1) = dmax;
    ne(k,1) = nec(points(ii));
    permind(k,1) = points(ii);
  end
  [dmin,ii] = min(ndc(points));
  if ~isempty(dmin)
    nd(k,2) = dmin;
    ne(k,2) = nec(points(ii));
    permind(k,2) = points(ii);
  end
  [dmed] = median(ndc(points));
  if ~isempty(dmed)
    nd(k,3) = dmed;
    [mm,ii] = min(abs(ndc(points)-dmed));
    ne(k,3) = nec(points(ii));  
    permind(k,3) = points(ii);
  end
  
end
%plot(nd,ne,'o',ndc,nec,'+');figure(1);