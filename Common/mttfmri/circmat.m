function y = circmat(x,n);

% Usage: y = circmat(x,n);
%  
%  x = vector
%  n = number of circular shifts.
%
%  y = output matrix where columns are n-circulant shifts of x.

x1 = x(:)';
N = length(x1);
if (n > N)
	error('circmat: n needs to be less than length of x');
end
    
x1 = [x1(1) x1(N:-1:2)];
y = toeplitz(x,x1(:,1:n));
