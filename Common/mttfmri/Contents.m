% 
% MATLAB Toolbox for Multiple Trial Type Event-Related fMRI  (mtfmri)
%
% Version 1.0
% Release Date: March 26, 2004
% Send Comments and Questions to ttliu@ucsd.edu
% See http://fmriserver.ucsd.edu/ttliu for updates.
% Copyright (c) 2004 UC Regents
%
%  Introduction
%  ------------
%  This toolbox is a companion to the papers
%  Liu and Frank, NeuroImage 21:387-400, 2004
%  Liu, NeuroImage, 21:401-413, 2004
%  Liu , "Errata for NeuroImage 21:387-400, 2004" coming soon 
%       (in the meantime, PDF file of errata is available at the website).
%
%  Overview
%  --------
%  m-sequences:           Use the gen_mseq function (see multsim3.m for examples)
%  clustered m-sequences: Use the gen_mseq function to make an
%                         m-sequence. Then use the clumpvec function
%                         to cluster the sequence. See multsim3.m and
%                         make_cluster.m for examples
%
%  permuted block:        Use gen_block and permute_block 
%                         See also multsim3 and make_mixed (these are
%                         older functions that don't use gen_block
%                         and permute_block, which are newer functions created from
%                         code in multsim3)
%
%  mixed design:          See make_mixed.m for example
%
%  *** To evaluate the performance of a design **** 
%  Use calc_meffdet to calculate efficiencies and detection powers. 
%  Use calcent or calcentvec to calculate conditional entropies.
%  Use tcurve to make up theoretical curves for efficiency and power. 
%
%  Getting Started
%  ---------------
%  Try out the scripts:
%
%    examp_blockperm -- example of block permutation
%    examp_mcluster  -- example of m-sequence clustering
%    meffdet_examp   -- an older example script with some useful
%                       printing functions
%
%  Example Scripts for the Paper
%  -----------------------------
%  1. Use a batch*.m file to generate some results
%  2. Use a plot*.m file to plot out the results
%  Example: batch_basis
%           plot_basis_effect
%
%  To see what files need to be executed prior to each plot*.m file,
%  look at the help info for that file. 
%  Example: batch01, make_mixed, and plot_nipaper01_ent need to be run
%           prior to plot_nipaper01.m
%
%  OPTION: You can also download example *.mat files from the website
%          and just plot out some pre-computed sample results. 
%  NOTE: Although the example scripts correspond to the 
%        figures in the papers, in some cases the number of 
%        permutations or paths has been decreased to reduce the execution
%        time. Thus, the results will tend to differ slightly from
%        the published results. In addition, the behavior of the random
%        number generator may differ from system to system. 
%
%
%  batch01.m        calls multsim3 to generate results for plot_nipaper01
%  batch_basis      calls multsim3 to generate results for plot_ms1b_basis2
%  batch_color      calls multsim3 to generate results for plot_ms1b_nc2
%  batch_nuisance   calls multsim3 to generate results for plot_nuisance
%  batch_cluster    example of using make_cluster
%
%  make_cluster     make up clustered and permuted block designs
%  make_mixed       make up example mixed designs

%
%
%  Figures from the Papers
%  -----------------------
%  plot_nipaper01     Figures  1 and 2 in Liu&Frank 2004
%  plot_nipaper01_ent used to generate entropy results for plot_nipaper_01
%                     (NOTE: this program can take a while to compute entropies).
%
%  plot_basis_effect  Figure 3 in Liu&Frank 2004
%  plot_colornoise    Figure 4 in Liu&Frank 2004
%
%  plot_blockpermute  Fig 1 in Liu 2004
%  plot_nuisance2     Figure 2  in Liu 2004
%  plot_cluster       Figures 3,4, 7 and 8 in Liu 2004
%  show_mlimit_ent    example of computing efficiencies and entropies
%                     of m-sequences (Fig 5 in Liu 2004)
%  plot_clusterexamp  show clustering example (Fig 6 in Liu 2004)
%
%
%  Design Programs
%  ---------------
%   multsim3         main file for running simulations 
%   calc_meffdet     calculate efficiencies and powers
%   gen_block        make up a block design
%   permute_block    permute a block design
%   clumpvec         algorithm for clumping a random sequence
%   optprob          returns optimum frequency of occurrence
%   tcurve           theoretical trade-off curves
%
%   calcent          calculate entropy of a sequence
%   calcent2         alternate way of calculating entropy
%   calcentvec       calculate entropies for matrix of sequences
%
%
%  M-sequence Programs
%  -------------------
%   gen_mseq         generate matrix of m-sequences, calls mseq2.m
%   mseq2            generate m-sequences (based on code from G. Buracas)
%   qadd             utility function for mseq2.m
%   qmult            utility function for mseq2.m
%   return_mtaps     returns number of possible m-sequences available. 
%
%  Utility programs
%  -----------------
%   approxtrace      approximation of the trace of the design matrix
%   arrow            arrow utility from E.A. Johhnson, included here for convenience
%   calc_freq        calculate frequencies of occurrence of trials in a  design
%   check_mdesign    load ASCII file and calculate efficiency and power
%   compress_clump   compress design permutations for display
%   dcon             utility function for contrasts
%   hemoresp         generate a gamma density function
%   laguerre         generate Laguerre basis functions
%   legendremat      generate matrix of legendre polynomials
%   printdesign      print out ASCII file of a design
%   shiftvec         shifts a vector
%   stimpatch        used to draw designs
%   voltbasis        generate gamma density basis functions

