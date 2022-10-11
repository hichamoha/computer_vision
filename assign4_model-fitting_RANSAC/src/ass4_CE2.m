% ------------ Computer Vision - Assignment 4 -------------------------
% ------------ Hicham Mohamad, hsmo@kth.se    -------------------------
%{
============   3 - Robust Homography Estimation and Stitching   =======

Computer Exercise 2. 
In this exercise you will use RANSAC to estimate homographies for 
creating panoramas.

You can use the two images a.jpg and b.jpg (see Figure 2). 
You will need to use VLfeat as in assignment 2 to generate potential 
matches, and then determine inliers using RANSAC. 

Begin by loading the two images in Matlab and displaying them. 
The images are partly overlapping. The goal is to place them on top 
of each other as in Figure 2. Use VLFeat to compute SIFT features 
for both images and match them.

How many SIFT features did you ?nd for the two images, respectively? 
How many matches did you ?nd? 
%}
%% load the data
clear all, close all
A = imread('a.jpg');
B = imread('b.jpg');

%% Use VLFeat to compute SIFT features for both images and match them
%{
The vector fA contains 4 rows. The first two rows are the coordinates of 
the detected features. The second two contains an orientation and scale 
for which the the feature was detected.
%}
% compute the features for image A and image B
[fA dA] = vl_sift( single(rgb2gray(A)) ); 
[fB dB] = vl_sift( single(rgb2gray(B)) );

% match the descriptors 
%{
Thus for each descriptor in image 1 the function finds the two best 
matches. If the quality of the best and thesecond best match is similar 
then the match is rejected, otherwise the best match is selected.
%}
matches = vl_ubcmatch(dA,dB);

% We extract matching points using:
xA = fA(1:2,matches(1,:)); 
xA3 = [xA;ones(1,size(xA,2))];
xB = fB(1:2,matches(2,:));
xB3 = [xB;ones(1,size(xA,2))];
[mA,nA] = size(xA);
[mA3,nA3] = size(xA3);

%% Plot the features together with the images
figure
imagesc(A)
hold on
vl_plotframe(fA);
title('Features together with image A')

% Compute the features for the second image 
%[fB dB] = vl_sift( single(rgb2gray(imCube2)), 'PeakThresh', 1);
figure
imagesc(B)
hold on
vl_plotframe(fB);
title('Features together with image B')

%% The following code randomly selects 10 matches, plots the two images 
% next to each other and plots lines between the matching points.
perm = randperm(size(matches ,2));
figure;
imagesc ([A B]);
hold on;
plot([xA(1,perm (1:10));  xB(1,perm (1:10))+ size(A ,2)], ...
     [xA(2,perm (1:10));  xB(2,perm (1:10))] ,'-');
hold  off;
title('Lines between the 10 matching points')
% How many of the matches appear to be correct?

%% Finding the homography H between the two images
%{
Now you should find a homography describing the transformation between 
the two images. Because not all matches are correct, you need to 
use RANSAC to find a set of good correspondences (inliers). 
To estimate the homography use DLT with a minimal number of points 
needed to estimate the homography. (Note that in this case the least 
squares system will have an exact solution, so normalization does 
not make any difference.) A reasonable threshold for inliers is 5 pixels.
 How many inliers did you find?

Next transform the images to a common coordinate system using the 
estimated homography.
%}

% Estimating homographies using RANSAC

iters = 10;
H = {};
inliersNbr = zeros(iters,1);
inliersIndx = {iters; 1};
for i=1:iters
    % pseudo random integers
    randind = randi(nA3, [1 4]);
    
    % Set up the DLT equations for matrix M
    % Note Xi is a 4x1 vector
    % each 0 on the left hand side represents a 1x4 block of zeros
    x1 = xA3(:,randind);
    x2 = xB3(:,randind);
    [m,n] = size(x1);
    
    M = zeros(n*3, m*3 + n);
    for iM=1:n
        M(3*(iM-1)+1,1:3) = x1(:,iM)';
        M(3*(iM-1)+2,4:6) = x1(:,iM)';
        M(3*(iM-1)+3,7:9) = x1(:,iM)';
        M(3*(iM-1)+1:3*iM,9+iM) = -x2(:,iM);
    end
    
    %Computes the singular value decomposition of M
    [U,S,V] = svd(M);
    
    % to find an eigenvector corresponding to the smallest eigenvalue 
    % we should select the last column of V
    vstar = V(:,end);
    H{i} = reshape(vstar(1:9),[3 3])';

    indexSet = [];
    %for j=1:n   %XXXXXXX correction
    for j=1:length(xA)    
        % Finds the the indices for which the distance to 
        % the plane is less than 0.1. 
        % Note: Works only if the 4th coordinate of all the 
        % points in X is 1.
        if(norm(pflat(H{i}*xA3(:,j)) - xB3(:,j)) <= 5)
            inliersNbr(i) = inliersNbr(i) + 1 ; 
            indexSet = [indexSet j];
        end
    end
    inliersIndx{i} = indexSet;
    
end 

%%  select the solution that gives the largest consensus set
largestConsensusIdx = find(inliersNbr == max(inliersNbr));
consensusSize = inliersNbr(largestConsensusIdx);
%disp(['Number of inliers with RANSAC ',num2str(consensusSize)])

bestH = H{largestConsensusIdx}

%% Image Stitching using the estimated H
% Next transform the images to a common coordinate system 
% using the estimated homography H.

tform = maketform('projective',bestH'); 
%Creates a transfomation that matlab can use for images 
%Note: imtransform uses the transposed homography 
transfbounds = findbounds(tform ,[1 1; size(A,2) size(A,1)]); 
%Finds the bounds of the transformed image 
xdata = [min([transfbounds(:,1); 1]) max([transfbounds(:,1); size(B,2)])]; 
ydata = [min([transfbounds(:,2); 1]) max([transfbounds(:,2); size(B,1)])]; 
%Computes bounds of a new image such that both the old ones will fit.

[newA] = imtransform(A,tform ,'xdata',xdata ,'ydata',ydata); 
%Transform the image using bestH

tform2 = maketform('projective',eye(3)); 
[newB] = imtransform(B,tform2 ,'xdata',xdata ,...
                     'ydata',ydata ,'size',size(newA)); 
%Creates a larger version of B
newAB = newB; 
newAB(newB < newA) = newA(newB < newA); 
%Writes both images in the new image. 
%(A somewhat hacky solution is needed 
%since pixels outside the valid image area are not allways zero...)

figure
imagesc(newAB)
title('CE2: The created panorama after estimating H')

%% Useful matlab commands:
%{
[fA dA] = vl_sift( single(rgb2gray(A)) ); 
[fB dB] = vl_sift( single(rgb2gray(B)) );

matches = vl_ubcmatch(dA,dB);

xA = fA(1:2,matches(1,:)); 
xB = fB(1:2,matches(2,:));

tform = maketform(?projective?,bestH ?); 
%Creates a transfomation that matlab can use for images 
%Note: imtransform uses the transposed homography 
transfbounds = findbounds(tform ,[1 1; size(A,2) size(A,1)]); 
%Finds the bounds of the transformed image 
xdata = [min([transfbounds(:,1); 1]) max([transfbounds(:,1); size(B,2)])]; 
ydata = [min([transfbounds(:,2); 1]) max([transfbounds(:,2); size(B,1)])]; 
%Computes bounds of a new image such that both the old ones will fit.

[newA] = imtransform(A,tform ,?xdata?,xdata ,?ydata?,ydata); 
%Transform the image using bestH

tform2 = maketform(?projective?,eye(3)); 
[newB] = imtransform(B,tform2 ,?xdata?,xdata ,?ydata?,ydata ,?size?,size(newA)); 
%Creates a larger version of B
newAB = newB; newAB(newB < newA) = newA(newB < newA); 
%Writes both images in the new image. 
%(A somewhat hacky solution is needed 
%since pixels outside the valid image area are not allways zero...)

%}
