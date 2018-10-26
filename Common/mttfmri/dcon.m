function dc = dcon(i,j,Q)
  dc = zeros(1,Q);
  if i == j
    dc(i) = 1;
  else
    dc(i) = 1;dc(j) = -1;
  end