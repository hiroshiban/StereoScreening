
To generate design files for multiple subjects called generateDesign with the following parameters:

 Subjects - number of subjects to generate designs for
 Designs - number of design files per subject
 Conditions - number of conditions in the experiment (including fixation)
 Reps - number of repeats of a sequence of conditions
 lookBack - number of previous trials which form a combination which is equally likely to appear in the trial sequence
 outputDir - directory to create design files into - this directory must exist

 Note: if the design file already exists then it will be overwritten if the same output directory is specified twice

 example:
 
 generateDesign(5,10,8,2,1,'d:\designs')
 
 will create 5 subject folders named 'subject 1' etc for each subject there will be 10 design files named design1.txt etc. in the relevant design folder.