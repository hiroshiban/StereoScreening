function y = rcircmat(x,n);
  
%   function y = rcircmat(x,n);
% 
%  generate circmat in the reverse direction of circmat
  
x1 = flipud(x(:));
x1 = circmat(x1,n);
y = flipud(x1);
