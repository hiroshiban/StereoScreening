function y = shiftvec(x,offset);


% function y = shiftvec(x,offset);
% shift a vector by an offset,
% positive = shift to the right and zero fill on the left.
% negative = shift to the left and zero fill on the right.

x = x(:);xlen = length(x);
if offset > 0
	y = [zeros(offset,1);x(1:(xlen-offset))];
else
	y = [x((1-offset):xlen);zeros(-offset,1)];
end
