function runpost(file,forw)

%runpost(file, forw)
%  where file is the basname of the files to be prprocessed.
%  forw is the forward model to be used.
% 

eval(['!post ' file ' ' forw]);
inputfil = file ;
plotsaga