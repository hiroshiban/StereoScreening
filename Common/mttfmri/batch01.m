% example batch file
% for showing the performance of designs
% similar to that in Fig 1 of Liu&Frank 2004

npts = 240;
numperm = 200;
numpermrand = 100; % in the paper, this was 1000, but here
                   % we make it 100 to keep simulation time down.
docolornoise = 0;
dobasis = 0;
doclump = 0;
for nevents = 2:5;
  multsim3;
end

