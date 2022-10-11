% ------------ Computer Vision - Assignment 2 -------------------------
% ------------ Hicham Mohamad, hsmo@kth.se
%{
2 - Calibrated vs Uncalibrated Reconstruction =========================
 ======================================================================

Computer Exercise 1

The file compEx1data.mat contains the 3D points of the 
reconstruction \textbf{X}, the camera matrices $P$, the image points x 
and the filenames \texttt{imfiles} of the images. 

Here \textbf{X} is a $4 \times 9471$ matrix containing the homogeneous 
coordinates for all 3D points, $x{i}$ is a $3 \times 9471$ matrix 
containing the homogeneous coordinates of the image points seen in image 
$i$ (NaN means that the point has not been detected in this image). 
$P{i}$ contains the camera matrix of image $i$ and $imfiles{i}$ contains 
the name of that image.
%}
clear all, close all
load('compEx1data.mat')

%% Displays the image
im = imread(imfiles{1});
figure
imagesc(im)    
colormap gray

%% Plot the 3D points of the reconstruction. 
figure
plot3(X(1,:),X(2,:),X(3,:),'.','Markersize',2);

title('CE1: 3D Points of the reconstruction')
xlabel('x'), ylabel('y'), zlabel('z')

%axis equal % axes have the same scale.
hold on

% Using the file plotcams.m to plot the cameras in the same figure. 
plotcams(P)
% Does this look like a reasonable reconstruction?  
% use axis equal otherwise you may get additionaldistortion
legend('3D points','Cameras')
axis equal

%% Project the 3D points into one of the cameras 
% only those that have been detected in that camera
i = 1;
Pi = P{i};
%xi = x{i};

xi = Pi*X;

% divide the points xi by the third coordinate before plotting 
[m,n] = size(xi);
xproj = zeros(m,n);
for j = 1 : n
    xproj(:,j) = pflat(xi(:,j)); 
end

%% Plot the image, the projected points, 
% and the image points in the same figure.
im = imread(imfiles{i});
figure
imagesc(im)    % Displays the image
colormap gray
hold on

% Determines which of the points are visible in image 
visible = isfinite(x{i}(1 ,:));

%Plots a '*' at each  point  coordinate
plot(x{i}(1, visible), x{i}(2, visible),'*');

% Plots a red 'o' at each  visible  point in xproj
plot(xproj(1,visible), xproj(2,visible),'ro');

title('Projected points and image points')
legend('Image points', 'Projected points')

%% New projective solutions using the two projective transformations
T1 = [1 0 0 0; 0 4 0 0; 0 0 1 0; 1/10 1/10 0 1];
T2 = [1 0 0 0; 0 1 0 0; 0 0 1 0; 1/16 1/16 0 1];

% modify the 3D points and cameras
% 3D points
Xtilde1 = T1*X;
Xtilde2 = T2*X;
% divide the points by the fourth coordinate before plotting
[m,n] = size(Xtilde1);
homXtilde1 = zeros(m,n);
homXtilde2 = zeros(m,n);
for j = 1 : n
    homXtilde1(:,j) = pflat(Xtilde1(:,j));
    homXtilde2(:,j) = pflat(Xtilde2(:,j));
end

%% modify cameras
Ptilde1 = {};
Ptilde2 = {};
for j=1:9
    Ptilde1{j} = P{j}/T1; % P * inverse of T1
    Ptilde2{j} = P{j}/T2;
end

%% Plot the 3D points and cameras in the same figure 
% for new solutions using T1
figure
plot3(homXtilde1(1,:),homXtilde1(2,:),homXtilde1(3,:),'.','Markersize',2);
hold on

% Using the file plotcams.m to plot the cameras in the same figure. 
plotcams(Ptilde1)

axis equal

title('CE1: Modified 3D Points using T1')
xlabel('x'), ylabel('y'), zlabel('z')
%axis equal % axes have the same scale.

%% for new solutions using T2
figure
plot3(homXtilde2(1,:),homXtilde2(2,:),homXtilde2(3,:),'.','Markersize',2);
hold on

% Using the file plotcams.m to plot the cameras in the same figure. 
plotcams(Ptilde2)

axis equal

title('CE1: Modified 3D Points using T2')
xlabel('x'), ylabel('y'), zlabel('z')

%% Project the new 3D points into one of the cameras
% only  those that have been detected in that camera
i = 1;
PT1i = Ptilde1{i};
PT2i = Ptilde2{i};
%xi = x{i};

xT1i = PT1i*Xtilde1;
xT2i = PT2i*Xtilde2;

% divide the points xi by the third coordinate before plotting 
[m1,n1] = size(xT1i);
xT1proj = zeros(m1,n1);
xT2proj = zeros(m1,n1);
for j = 1 : n1
    xT1proj(:,j) = pflat(xT1i(:,j)); 
    xT2proj(:,j) = pflat(xT2i(:,j));
end

%% Plot the image, the projected points, using T1
% and the image points in the same figure.
im1 = imread(imfiles{i});
figure
imagesc(im1)    % Displays the image
colormap gray
hold on

% Determines which of the points are visible in image 
%visible = isfinite(x{i}(1 ,:));

%Plots a '*' at each  point  coordinate
plot(x{i}(1, visible), x{i}(2, visible),'*');

% Plots a red 'o' at each  visible  point in xproj
plot(xT1proj(1,visible), xT1proj(2,visible),'ro');

title('New Projected points and image points using T1')
legend('Image points', 'Projected points')

%% Plot the image, the projected points, using T2
% and the image points in the same figure.
im1 = imread(imfiles{i});
figure
imagesc(im1)    % Displays the image
colormap gray
hold on

% Determines which of the points are visible in image 
%visible = isfinite(x{i}(1 ,:));

%Plots a '*' at each  point  coordinate
plot(x{i}(1, visible), x{i}(2, visible),'*');

% Plots a red 'o' at each  visible  point in xproj
plot(xT2proj(1,visible), xT2proj(2,visible),'ro');

title('New Projected points and image points using T2')
legend('Image points', 'Projected points')

%% Exercise 3
K3 = [320 0 320; 0 320 240; 0 0 1];
A = [0; 240; 1]; B = [640; 240; 1];
Atilde = K3\A
Btilde = K3\B

%% Exercise 4
% consider the camera
P4 = [1000 -250 250*sqrt(3) 500; 
      0 500*(sqrt(3)-(1/2)) 500*(1+(sqrt(3)/2)) 500;
      0 -1/2 sqrt(3)/2 1]; 
K = [1000 0 500; 0 1000 500; 0 0 1]; 
calP4 = K\P4; 

% Normalize the corners of images of size 1000x1000 
leftup = [0;0;1];
rightup = [0;1000;1];
leftdown = [1000;0;1];
rightdown = [1000;1000;1];
center = [500; 500; 1]

leftupTilde = K\leftup;
rightupTilde = K\rightup;
leftdownTilde = K\leftdown;
rightdownTilde = K\rightdown;
centerTilde = K\center;

%% 4 - Q Factorization and Computation of K ==========================
% Exercise 5
d = sqrt(0 + (1400)^2 + 0)

%% Computer Exercise 2
% ====================================================================
[K1,R1] = rq(PT1i);
[K2,R2] = rq(PT2i);

%% Make sure that element K(3,3) = 1 by division;
K1test = K1./K1(3,3)
K2test = K2./K2(3,3)
% Do they represent the same transformation?  
% Two matrices can differ by a scale factor and still give 
% the same transformation.

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