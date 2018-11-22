function nn = members(dat)
% individuals = members(dat)
% returns all individual members of dat in sorted order as a vector

nn=sort(dat);
l = length(nn);
i = 1;
while i<l
  s=nn(i);
  f = find(nn(:) == s);
  f(find(f(:) == i)) = [];
  nn(f) = [];
  i=i+1;
  l=length(nn);
end

