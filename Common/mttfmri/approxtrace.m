function M = approxtrace(k,m,N);

x = 0:(k-1);
M = sum((1 - (1-x/N)*m/N).*(1-x/N)*m);
