function [h,dh] = hemoresp(t,tau,n,delay);

% MODS
%  020516 TTL  added temporal derivative as an output

tdel = t-delay;
h = zeros(size(t));
dh = zeros(size(t));
pos_t = (find(tdel >=0));
tdel = tdel(pos_t);

h(pos_t) = (tdel/tau).^n.*exp(-tdel/tau)/tau/prod(1:n);
dh(pos_t) = n*tdel.^(n-1)/tau^n.*exp(-tdel/tau)/tau/prod(1:n) - ...
            1/tau*(tdel/tau).^n.*exp(-tdel/tau)/tau/prod(1:n);
h = h(:);dh = dh(:);
