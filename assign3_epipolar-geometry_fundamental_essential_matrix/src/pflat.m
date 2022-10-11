%{
Computer Exercise1.
pflat() divides the homogeneous coordinates with their last entry for 
points of any dimensionality. (we assume that none of the points have 
last homogeneous coordinate zero.)
%}

function x = pflat(hmgCoord)
x = hmgCoord ./ hmgCoord(end,:); 
end

