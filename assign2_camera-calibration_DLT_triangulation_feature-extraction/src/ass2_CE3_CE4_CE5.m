% ====================================================================
%      Assignment 2: Computer Exercises 3-4-5
%      Hicham Mohamad hsmo@kth.se
% ====================================================================

% ===========  5 - Direct Linear Transformation DLT ==================
%{
Figure 2 shows two imagescube1.jpg and cube2.jpgof a scene with a 
Rubics cube. The file compEx3data.matcontains a point model Xmodel of 
the visible cube sides, the measured projections x of the model points 
in the two images and two variables startind, endind that can be used 
for plotting lines on themodel surface.
%}
clear all, close all
load('compEx3data.mat')

%% Normalize the measured points by applying a transformation N  
% transformation N that subtracts the mean of the points and then
% re-scales the coordinates by the standard deviation in each coordinate.

% mean and std of the points in x1
% Computes the mean value of the 1st and 2nd rows of x{i}
meanx1 = mean(x{1}(1:2,:),2);
%Standard deviation of the 1st and 2nd rows of x{i}
stdx1 = std(x{1}(1:2,:),0,2);

% mean and std of the points in x2
meanx2 = mean(x{2}(1:2,:),2);
stdx2 = std(x{2}(1:2,:),0,2);

% transformation N
N1 = [1/stdx1(1) 0 -1/stdx1(1)*meanx1(1);
      0 1/stdx1(2) -1/stdx1(2)*meanx1(2);
      0 0 1];
N2 = [1/stdx2(1) 0 -1/stdx2(1)*meanx2(1);
      0 1/stdx2(2) -1/stdx2(2)*meanx2(2);
      0 0 1];
% Normalizing the measured points x
xtilde1 = N1*x{1};
xtilde2 = N2*x{2};
 
%% Plot the normalized points in a new figure
figure
plot(xtilde1(1,:), xtilde1(2,:),'.')
title('The normalized points with transformation N1')
axis equal

figure
plot(xtilde2(1,:), xtilde2(2,:),'.')
title('The normalized points with transformation N2')
axis equal
% Does it look like the points are centered around(0,0)with mean 
% distance 1 to(0,0)?

%% Set up the DLT equations for resectioning: M1 and M2
% Note Xi is a 4x1 vector
% each 0 on the left hand side represents a 1x4 block of zeros
Xnb = size(Xmodel,2);

M1 = zeros(Xnb*3, 4*3+Xnb);
for iM=1:Xnb
    M1(3*(iM-1)+1,1:4) = [Xmodel(:,iM); 1]';
    M1(3*(iM-1)+2,5:8) = [Xmodel(:,iM); 1]';
    M1(3*(iM-1)+3,9:12) = [Xmodel(:,iM); 1]';
    M1(3*(iM-1)+1:3*iM,12+iM) = -xtilde1(:,iM);
end

M2 = zeros(Xnb*3, 4*3+Xnb);
for iM=1:Xnb
    M2(3*(iM-1)+1,1:4) = [Xmodel(:,iM); 1]';
    M2(3*(iM-1)+2,5:8) = [Xmodel(:,iM); 1]';
    M2(3*(iM-1)+3,9:12) = [Xmodel(:,iM); 1]';
    M2(3*(iM-1)+1:3*iM,12+iM) = -xtilde2(:,iM);
end

%% solve the resulting homogeneous least squares system using SVD
%Computes the singular value decomposition of M
[U1,S1,V1] = svd(M1);
[U2,S2,V2] = svd(M2);

% S'S is a diagonal matrix and contains the eigenvalues
% Is the smallest eigenvalue close to zero? 
diag1 = S1'*S1;
diag2 = S2'*S2;

% to find an eigenvector corresponding to the smallest eigenvalue 
% we should select the last column of V
v1star = V1(:,end);
v2star = V2(:,end);

% How about||Mv||?
M1v1 = norm(M1*v1star);
M2v2 = norm(M2*v2star);

%% Extract the entries of the camera from the solution 
% Make sure that you select the solution where the points are 
% in front of the camera. (If X has 4th coordinate 1 then the 
% 3rd coordinate of PX should be positive for X to be in front 
% of the camera.)

P1tilde = reshape(-v1star(1:12),[4 3])';
% set up the camera matrix
P1 = N1\P1tilde;
% Project the model points into the images
P1Xmodel = P1*[Xmodel;ones(1,size(Xmodel,2))];

P2tilde = reshape(-v2star(1:12),[4 3])';
% set up the camera matrix
P2 = N2\P2tilde;
% Project the model points into the images
P2Xmodel = P2*[Xmodel;ones(1,size(Xmodel,2))];

%% Project the model points into the images
% We first transform the camera matrix to the original 
% (unnormalized) coordinate system, as in Exercise 7, before projecting.
% divide the points xi by the third coordinate before plotting 
[m3,n3] = size(P1Xmodel);
xP1proj = zeros(m3,n3);
xP2proj = zeros(m3,n3);
for j = 1 : n3
    xP1proj(:,j) = pflat(P1Xmodel(:,j)); 
    xP2proj(:,j) = pflat(P2Xmodel(:,j));
end

%% load the images cube1 and cube2
imCube1 = imread('cube1.jpg');
imCube2 = imread('cube2.jpg');

%% Plot the measured image points in the same figure
figure;
imagesc(imCube1)
hold on
plot(x{1}(1,:),x{1}(2,:),'*')
hold on;
plot(xP1proj(1,:), xP1proj(2,:),'ro');

title('Projected points with P1 and image points')
legend('Image points', 'Projected points')

figure;
imagesc(imCube2)
hold on
plot(x{2}(1,:),x{2}(2,:),'*')
hold on;
plot(xP2proj(1,:), xP2proj(2,:),'ro');

title('Projected points with P2 and image points')
legend('Image points', 'Projected points')
% Are they close to each other?

%% plot the camera centers and viewing directions in the same plot 
% as the model points. Does the result look reasonable? 
figure;
plot3(Xmodel(1,:),Xmodel(2,:),Xmodel(3,:),'.','Markersize',6)
hold on;
%Plots the lines of the cube model 
%(works  only if all points are included)
plot3([ Xmodel(1,startind );  Xmodel(1,endind )],...
      [Xmodel(2,startind );  Xmodel(2,endind )],...
      [Xmodel(3,startind );  Xmodel(3,endind)],'r-');
axis equal

plotcams({P1, P2})

title('CE3: 3D-points and camera centers')
xlabel('x'), ylabel('y'), zlabel('z')

%% Compute the inner parameters of the first camera using rq.m 
[K1cube,R1cube] = rq(P1);
% How can we know that these are the "true" param-eters? 
% Why is there no ambiguity as in Exercise 1?

%% 6  Feature Extraction and Matching using SIFT =======================
% ======================================================================
%{
First we download and start VLFeat. 
Go tohttp://www.vlfeat.org/download.htmland 
extractthe binary package to a directory of your choice. 
Then start Matlab, go the
H:\vlfeat\toolbox subdirectory and run vlsetup. 
Now you should see the following message:
** Welcome to the VLFeat Toolbox **
We will now be able to use VLFeat throughout this Matlab session.
First load the imagescube1.jpgandcube2.jpgfrom Exercise 3.
Compute sift features using vlfeat.
%}
close all
% run vlsetup
%C:\Users\hashu\Downloads\vlfeat-0.9.21\toolbox\vl_setup.m

% First load the imagescube1.jpgandcube2.jpg
imCube1 = imread('cube1.jpg');
imCube2 = imread('cube2.jpg');

%% Compute sift features using vlfeat
%{
The SIFT detector searches for peaks in scale space (similar to peaks in 
the autocorrelation function, see lecture notes). The second argument 
filters out peaks that are too small.
The vector f1 contains 4 rows. The first two rows are the coordinates of 
the detected features. The second two contains an orientation and scale 
for which the the feature was detected.
%}
[f1 d1] = vl_sift( single(rgb2gray(imCube1)), 'PeakThresh', 1);

%% Plot the features together with the images
figure
imagesc(imCube1)
hold on
vl_plotframe(f1);
title('Features together with image 1')

% Compute the features for the second image 
[f2 d2] = vl_sift( single(rgb2gray(imCube2)), 'PeakThresh', 1);
figure
imagesc(imCube2)
hold on
vl_plotframe(f2);
title('Features together with image 2')

%% match the descriptors 
[matches ,scores] = vl_ubcmatch(d1,d2);

% Thus for each descriptor in image 1 the function finds the two best 
%matches. If the quality of the best and thesecond best match is similar 
%then the match is rejected, otherwise the best match is selected.

% We can now extract matching points using:
x1 = [f1(1,matches (1 ,:));f1(2,matches (1 ,:))];
x2 = [f2(1,matches (2 ,:));f2(2,matches (2 ,:))];

%% The following code randomly selects 10 matches, plots the two images 
% next to each other and plots lines between the matching points.
perm = randperm(size(matches ,2));
figure;
imagesc ([imCube1 imCube2]);
hold on;
plot([x1(1,perm (1:10));  x2(1,perm (1:10))+ size(imCube1 ,2)], ...
     [x1(2,perm (1:10));  x2(2,perm (1:10))] ,'-');
hold  off;
title('Lines between the 10 matching points')
% How many of the matches appear to be correct?

%% =========== 7  Triangulation using DLT ===========================
%====================================================================
% Using the estimated cameras from Computer Exercise 3 you will now 
% triangulate the points detected in Computer Exercise 4.

% Set up the DLT equations for triangulation
X = [];
for j=1:size(x1,2)
    M = [P1 -[x1(:,j);1] [0 0 0]';
         P2 [0 0 0]' -[x2(:,j);1]];
     
    % and solve the homogeneous least squares system 
    % do this in a loop, once for each point
    [U,S,V] = svd(M);
    vstar = V(:,end);
    X = [X vstar(1:4,:)];
end

%% Project the computed points into the two images
P1X = P1*X;
P2X = P2*X;

[m4,n4] = size(P1X);
x1proj = zeros(m4,n4);
x2proj = zeros(m4,n4);
for j = 1 : n4
    x1proj(:,j) = pflat(P1X(:,j)); 
    x2proj(:,j) = pflat(P2X(:,j));
end

%% compare projected Triangulation points with the corresponding SIFT-points
% Plot the reconstructed image points in the same figure
figure;
imagesc(imCube1)
hold on
plot(x1(1,:),x1(2,:),'*')
hold on;
plot(x1proj(1,:), x1proj(2,:),'ro');

title('CE5: Projected Triangulation points and SIFT points into P1')
legend('SIFT-points', 'Triangulation points')

figure;
imagesc(imCube2)
hold on
plot(x2(1,:),x2(2,:),'*')
hold on;
plot(x2proj(1,:), x2proj(2,:),'ro');

title('CE5: Projected Triangulation and SIFT points into P2')
legend('SIFT-points', 'Triangulation points')

%% compare projected Triangulation points with the normalizaton points
% Plot the reconstructed image points in the same figure
figure;
imagesc(imCube1)
hold on
plot(xP1proj(1,:),xP1proj(2,:),'r*')
hold on;
plot(x1proj(1,:), x1proj(2,:),'go');

title('CE5: Projected Triangulation and normalization points into P1')
legend('Normalization points', 'Triangulation points')

figure;
imagesc(imCube2)
hold on
plot(xP2proj(1,:),xP2proj(2,:),'r*')
hold on;
plot(x2proj(1,:), x2proj(2,:),'go');

title('CE5: Projected Triangulation and normalization points into P2')
legend('Normalization points', 'Triangulation points')

%% Plot the remaining 3D points, the cameras and the cube model 
% Finds the points with reprojection error 
% less than 3 pixels in both images
good_points = (sqrt(sum((x1 - x1proj(1:2 ,:)).^2)) < 3 &...
               sqrt(sum((x2 - x2proj(1:2 ,:)).^2)) < 3);

% Removes points that are not good enough
goodX = X(:, good_points);
[m5,n5] = size(goodX);
homgoodX = zeros(m5,n5);
for j = 1 : n5
    homgoodX(:,j) = pflat(goodX(:,j)); 
end

%% Plot the reconstructed 3D points 
% and compare with the corresponding SIFT-points
figure;
plot3(homgoodX(1,:),homgoodX(2,:),homgoodX(3,:),'.','Markersize',2)
hold on;
xlabel('x'), ylabel('y'), zlabel('z')

% Using the file plotcams.m to plot the cameras in the same figure. 
plotcams({P1,P2})
hold on
%Plots the lines of the cube model 
%(works  only if all points are included)
plot3([Xmodel(1,startind); Xmodel(1,endind)],...
      [Xmodel(2,startind); Xmodel(2,endind)],...
      [Xmodel(3,startind); Xmodel(3,endind)],'-');
%grid on
axis equal
title('CE5: The reconstructed 3D-points and camera centers')

%% Useful Matlab commands
%{
im = imread(imfiles{i});
%Reads the imagefile with name in imfiles{i}

visible = isfinite(x{i}(1 ,:));
% Determines which of the points are visible in image 

iplot(x{i}(1, visible), x{i}(2, visible),'*');
%Plots a '*' at each  point  coordinate

plot(xproj(1,visible), xproj(2,visible),?ro?);
%Plots a red 'o' at each  visible  point in xproj

plot3(X(1,:),X(2,:),X(3,:),'.','Markersize' ,2);
%Plots a small  '.' at all the 3D points.

Useful  matlab  commands: CE3
%Computes the mean value of the 1st and 2nd rows of x{i}
mean(x{i}(1:2 ,:) ,2)

%Standard deviation  of the 1st and 2nd rows of x{i}
std(x{i}(1:2 ,:) ,0 ,2)

%Computes the singular value decomposition of M
[U,S,V] = svd(M);

%Takes the first 12 entries of sol and row -stacks them in a 3x4 matrix
P = reshape(sol (1:12) ,[4  3])';

%Plots the lines of the cube model 
(works  only if all points are included)
plot3([ Xmodel(1,startind );  Xmodel(1,endind )],...
      [Xmodel(2,startind );  Xmodel(2,endind )],...
      [Xmodel(3,startind );  Xmodel(3,endind)],'b-');

% Computes  the  projections
xproj1 = pflat(Pa*X);
xproj2 = pflat(Pb*X);

% Finds the points with reprojection error less than 3 pixels 
in both images
good_points = (sqrt(sum((x1-xproj1 (1:2 ,:)).^2))  < 3 &...
               sqrt(sum((x2-xproj2 (1:2 ,:)).^2))  < 3);

% Removes points that are not good enough
X = X(:, good_points );

%}