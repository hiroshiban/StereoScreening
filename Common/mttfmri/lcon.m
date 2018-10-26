
function lc = lcon(i,j,Q,k);
  lc = kron(dcon(i,j,Q),eye(k));
