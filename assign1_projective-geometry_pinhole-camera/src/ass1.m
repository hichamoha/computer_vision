% ------------ Computer Vision - Assignment 1 -------------------------
% ------------ Hicham Mohamad, hsmo@kth.se
%{
Points in Homogeneous Coordinates - Computer Exercise 1

Write a matlab function pflat that divides the homogeneous coordinates 
with their last entry for points of any dimensionality. (You may assume 
that none of the points have last homogeneous coordinate zero.)
Apply the function to the points in x2D and x3D in the file compEx1.mat, 
and plot the result.
%}
clear all
load('compEx1.mat')

% Apply the function to the points in x2D
[m2,n2] = size(x2D);
P2D = zeros(m2,n2);
for i = 1 : n2
    P2D(:,i) = pflat(x2D(:,i)); 
end

% Plot the result for the points in x2D
figure
plot(P2D(1,:),P2D(2,:),'.')
title('CE1: Points in x2D')
xlabel('x'), ylabel('y')
axis equal % axes have the same scale.

% Apply the function to the points in x3D
[m3,n3] = size(x3D);
P3D = zeros(m3,n3);
for i = 1 : n3
    P3D(:,i) = pflat(x3D(:,i)); 
end

% Plot the result for the points in x3D
figure
plot3(P3D(1,:),P3D(2,:),P3D(3,:),'.')
title('CE1: Points in x3D')
xlabel('x'), ylabel('y'), zlabel('z')
axis equal % axes have the same scale.

%% 3 - Lines - Computer Exercise 2
% In the file compEx2.mat there are three pairs of image points.  
im = imread('compEx2.JPG'); % Loads the image compEx2.JPG
load('compEx2')

figure
imagesc(im)    % Displays the image
colormap gray  % changes the colormap of the current image to gray scale
hold on        % Prevents the plot command from clearing the figure before 
               % plotting

% Plot the image points in the same figure as the image.
plot(p1(1,:), p1(2,:), 'or', ...
     p2(1,:), p2(2,:), 'xr', ...
     p3(1,:), p3(2,:), '+r', 'linewidth',1 )

% For each pair of points compute the line going through the points.
l1 = null(p1');
l2 = null(p2');
l3 = null(p3');

% Use the function rital to plot the lines in the same image.
rital(l1)
rital(l2)
rital(l3)
% Do these lines appear to be parallel (in 3D)?

%% Compute the point of intersection 
% between the second and third line (the lines obtained from the pairs p2 
% and p3).
x23 = null([l2,l3]')
cartesx23 = pflat(x23)

%Plot this point in the same image.
plot(cartesx23(1), cartesx23(2), '*w', 'linewidth', 1)

% Compute the distance between the first line and the intersection point.
% Is it close to zero? Why/why not?
d = (cartesx23' * l1) / norm(l1(1:2))

%% 4  Projective Transformations - Computer Exercise 3
% The filecompEx3.mat contains the start and end points of a set of lines.  
% Plotting the lines gives the grid in Figure 2
%clear all
load('compEx3')
[m,n] = size(startpoints);

%% Plots a blue line between each startpoint and endpoint
figure
plot([startpoints(1,:);  endpoints(1,:)], ...
     [startpoints(2,:);  endpoints(2,:)],'b-');

%% The projective mappings given by the matrices
H1 = [sqrt(3) -1 1; 1 sqrt(3) 1; 0 0 2];
H2 = [1 -1 1; 1 1 0; 0 0 1]; 
H3 = [1 1 0; 0 2 0; 0 0 1];
H4 = [sqrt(3) -1 1; 1 sqrt(3) 1; 1/4 1/2 2];

%% Transformations with projective mappings given by matrix H1
% compute the transformations of the given start and endpoints 
startTrans1 = H1 * [startpoints; ones(1, n)];
endTrans1 = H1 * [endpoints; ones(1,n)];

% compute cartesian coordinates by using pflat function
cartStartTrans1 = zeros(m+1,n);
cartEndTrans1 = zeros(m+1,n);
for i = 1 : n
    cartStartTrans1(:,i) = pflat(startTrans1(:,i)); 
    cartEndTrans1(:,i) = pflat(endTrans1(:,i));
end

%plot the lines between start and endpoints
figure
plot([cartStartTrans1(1,:);  cartEndTrans1(1,:)], ...
     [cartStartTrans1(2,:);  cartEndTrans1(2,:)],'b-');
 
% use the axis equal command, otherwisethe figures might look distorted
axis equal

%% Transformations with projective mappings given by matrix H2
% compute the transformations of the given start and endpoints 
startTrans2 = H2 * [startpoints; ones(1, n)];
endTrans2 = H2 * [endpoints; ones(1,n)];

% compute cartesian coordinates by using pflat function
cartStartTrans2 = zeros(m+1,n);
cartEndTrans2 = zeros(m+1,n);
for i = 1 : n
    cartStartTrans2(:,i) = pflat(startTrans2(:,i)); 
    cartEndTrans2(:,i) = pflat(endTrans2(:,i));
end

%plot the lines between start and endpoints
figure
plot([cartStartTrans2(1,:);  cartEndTrans2(1,:)], ...
     [cartStartTrans2(2,:);  cartEndTrans2(2,:)],'b-');
 
% use the axis equal command, otherwisethe figures might look distorted
axis equal

%% Transformations with projective mapping given by matrix H3
% compute the transformations of the given start and endpoints 
startTrans3 = H3 * [startpoints; ones(1, n)];
endTrans3 = H3 * [endpoints; ones(1,n)];

% compute cartesian coordinates by using pflat function
cartStartTrans3 = zeros(m+1,n);
cartEndTrans3 = zeros(m+1,n);
for i = 1 : n
    cartStartTrans3(:,i) = pflat(startTrans3(:,i)); 
    cartEndTrans3(:,i) = pflat(endTrans3(:,i));
end

%plot the lines between start and endpoints
figure
plot([cartStartTrans3(1,:);  cartEndTrans3(1,:)], ...
     [cartStartTrans3(2,:);  cartEndTrans3(2,:)],'b-');
 
% use the axis equal command, otherwisethe figures might look distorted
axis equal

%% Transformations with projective mappings given by matrix H4
% compute the transformations of the given start and endpoints 
startTrans4 = H4 * [startpoints; ones(1, n)];
endTrans4 = H4 * [endpoints; ones(1,n)];

% compute cartesian coordinates by using pflat function
cartStartTrans4 = zeros(m+1,n);
cartEndTrans4 = zeros(m+1,n);
for i = 1 : n
    cartStartTrans4(:,i) = pflat(startTrans4(:,i)); 
    cartEndTrans4(:,i) = pflat(endTrans4(:,i));
end

%plot the lines between start and endpoints
figure
plot([cartStartTrans4(1,:);  cartEndTrans4(1,:)], ...
     [cartStartTrans4(2,:);  cartEndTrans4(2,:)],'b-');
 
% use the axis equal command, otherwisethe figures might look distorted
axis equal
%{
Which of the transformations preserve lengths between  points?   
Which preserve angles between lines?   
Which maps parallel lines to parallel lines?  
Classify the transformations into euclidean, similarity, affine and 
projectivetransformations.
%}
%% 5 The Pinhole Camera 
% Exercise 5
P = [1 0 0 0; 0 1 0 0; 0 0 1 1];
% Compute the camera center (position) of the camera 
C = null(P) % computes the nullspace of P
c = pflat(C)

% computer the principal axis (viewing direction)
V = P(3 ,1:3) % extracts  elements P31 , P32 and P33

%% Computer Exercise 4
% Load and plot the images compEx4im1.jpg and compEx4im2.jpg
% The file compEx4.mat contains the camera matrices P1, P2 and a point 
% model U of the statue.
im1 = imread('compEx4im1.JPG'); % Loads the images
im2 = imread('compEx4im2.JPG');
load('compEx4')

%% Display the images im1 and im2
figure
imagesc(im1)    % Displays the image
colormap gray

figure
imagesc(im2)
colormap gray

%% Compute the camera centers and principal axes of the cameras
% Compute the camera center (position) of the camera 
C1 = null(P1) % computes the nullspace
c1 = pflat(C1)
C2 = null(P2)
c2 = pflat(C2)

% computer the principal axis (viewing direction)
V1 = P1(3 ,1:3) % extracts  elements P31 , P32 and P33
V2 = P2(3 ,1:3)

% normalize the principal axes
V1norm = V1./norm(V1)
V2norm = V2./norm(V2)

%% Plot the 3D-points in U and the camera centers in the same 3D plot 
% make sure that the 4th coordinate of U is one before you plotting
% Apply the function pflat() to the points in U
[mU,nU] = size(U);
Upoints = zeros(mU,nU);
for i = 1 : nU
    Upoints(:,i) = pflat(U(:,i)); 
end

% plot but with smaller points
figure
plot3(Upoints(1,:),Upoints(2,:),Upoints(3,:),'.','Markersize',2);
title('CE4: Points in U')
xlabel('x'), ylabel('y'), zlabel('z')
axis equal % axes have the same scale.
hold on

plot3(c1(1,:),c1(2,:),c1(3,:),'x','Markersize',10);
plot3(c2(1,:),c2(2,:),c2(3,:),'gx','Markersize',10);

% In addition plot a vector in the direction of the principal axes (
% viewing direction) from the camera center
quiver3(c1(1),c1(2),c1(3),V1norm(1),V1norm(2),V1norm(3),10)
% quiver3 Plots a vector v starting from
% the point a and rescales the sise by s
quiver3(c2(1),c2(2),c2(3),V2norm(1),V2norm(2),V2norm(3),10)

%% Project the points in U into the cameras P1 and P2 
U1 = P1*U;
U2 = P2*U;

[mU1,nU1] = size(U1);
u1 = zeros(mU1,nU1);
u2 = zeros(mU1,nU1);
for i = 1 : nU1
    u1(:,i) = pflat(U1(:,i)); 
    u2(:,i) = pflat(U2(:,i));
end

% plot the result in the same plots as the images
% Display the images im1 and im2
figure
imagesc(im1)    % Displays the image
colormap gray
hold on
plot(u1(1,:), u1(2,:),'.','Markersize',2)

figure
imagesc(im2)
colormap gray
hold on
plot(u2(1,:), u2(2,:),'.','Markersize',2)

% Does the result look reasonable?

%% useful  matlab  commands:
%repmat(a,[m n]) % Creates a block  matrix  with mn  copies  of a.
%hold on %Prevents the plot command from clearing the figure before plotting
%hold  off %Makes the plot command clear the figure before plotting
%null(A) %computes the nullspace of A



