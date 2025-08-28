% read the SHCs
function [nmax,cnm ,snm] = readSHC (file)
 n = file(: , 1);
 m = file(: , 2);
  for (i=1:size(file,1))
cnm(n(i)+1,m(i)+1) = file( i ,3 ) ;
snm(n(i)+1,m(i)+1) = file (i , 4 ) ;

  end      
nmax = file(end , 1);
end