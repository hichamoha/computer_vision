function rital(linjer,st)
%linjer is a 3xm matrix containing the homogeneous coordinates of n 2D
%lines
if nargin == 1,
 st='-';
end;
if size(linjer)==0,
  slask=[];
else
  [slask,nn]=size(linjer);
  rikt=psphere([linjer(2,:);-linjer(1,:);zeros(1,nn)]);
  punkter=pflat(cross(rikt,linjer));
  for i=1:nn;
   plot([punkter(1,i)-2000*rikt(1,i) punkter(1,i)+2000*rikt(1,i)], ...
        [punkter(2,i)-2000*rikt(2,i) punkter(2,i)+2000*rikt(2,i)],st);
  end;
  slask=[];
end;
