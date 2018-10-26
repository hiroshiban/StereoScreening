function pb = permute_block(basevec);

% Permute a block design
%  
%   pb = permute_block(basevec);
%
  
  npts = length(basevec);
  looping = 1;
  while looping % loop until we find two points with different values
    p1 = unidrnd(npts);
    p2 = unidrnd(npts);
    if basevec(p1) ~= basevec(p2)
      tmp = basevec(p1);
      basevec(p1)=basevec(p2);
      basevec(p2)=tmp;
      looping = 0;
    end
  end
  pb = basevec;