% ------------ Computer Vision - Assignment 3 -------------------------
% ------------ Hicham Mohamad, hsmo@kth.se
%{
2 - The Fundamental Matrix

In this exercise you will compute the fundamental matrix for the two 
images in Figure 1 of a part of the fort Kronan in Gothenburg. 
The file compEx1data.mat contains a cell x with matched points 
for the two images.
compute normalization matrices N_1 and N_2. These matrices should subtract the mean and re-scale using the standard 
% deviation, as in assignment 2.  Normalize the image points 
% of the two images with N_1 and N_2 respectively.
%}

% load the data
clear all, close all
load('compEx1data');
image1 = imread('kronan1.jpg');
image2 = imread('kronan2.jpg');

%% mean and std of the points in x1 and x2
% Computes the mean value of the 1st and 2nd rows of x{i}
meanx1 = mean(x{1}(1:2,:),2);
%Standard deviation of the 1st and 2nd rows of x{i}
stdx1 = std(x{1}(1:2,:),0,2);

% mean and std of the points in x2
meanx2 = mean(x{2}(1:2,:),2);
stdx2 = std(x{2}(1:2,:),0,2);

%% transformation matrices N1 and N2
% These matrices should subtract the mean and re-scale using the standard 
% deviation, as in assignment 2. 
N1 = [1/stdx1(1) 0 -1/stdx1(1)*meanx1(1);
      0 1/stdx1(2) -1/stdx1(2)*meanx1(2);
      0 0 1];
N2 = [1/stdx2(1) 0 -1/stdx2(1)*meanx2(1);
      0 1/stdx2(2) -1/stdx2(2)*meanx2(2);
      0 0 1];

%% Normalize the image points 
% of the two images with N_1 and N_2 respectively.
xn{1,1} = N1*x{1};
xn{2,1} = N2*x{2};

%% Set up the matrix M in the eight point algorithm 
% (use all the points), 
% Computes a 3x3 matrix containing all multiplications
% of coordinates from x1n(:,i) and x2n(:,i)
for i=1:2008
    xx = xn{2}(:,i)*xn{1}(:,i)';
    % Reshapes the matrix above and adds to the M matrix
    M(i,:) = xx(:)';
end

%% solve the homogeneous least squares system using SVD
[~,S,V] = svd(M);

% to find an eigenvector corresponding to the smallest eigenvalue 
% we should select the last column of V
vstar = V(:,end);

% Check that the minimum singular value and ||Mv|| are both small
% S'S is a diagonal matrix and contains the eigenvalues
% Is the smallest eigenvalue close to zero? 
D = S'*S;
% How about||Mv||?
Mv = norm(M*vstar);

%% Forms an F-matrix from the solution v of the least squares problem
Fn = reshape(vstar,[3  3]);
test_detFn = det(Fn);

% Resulting F may not have det(F) = 0. 
% Pick the closest matrix A with det(A) = 0.
% Can be solved using svd, in matlab:
[FU,FS,FV] = svd(Fn); 
FS(3,3) = 0;
% The resulting normalized fundamental matrix
A = FU * FS * FV';
% check the det(A)
test_detA = det(A);

%% Computes and plots all the epipolar constraints (should be roughly 0)
figure
plot(diag(xn{2}' * A * xn{1}),'.');
title('CE1: Plot of all the epipolar constraints')

%% Compute the un-normalized fundamental matrix F 
%{
Compute the un-normalized fundamental matrix F (using the formula from 
exercise 3) and the epipolar lines l=F x_1. 
Pick 20 points in the second image at random and plot these in the same 
figure as the image. 
Also plot the corresponding epipolar lines in the same image using the 
function rital.m. Are they close to each other?
%}

% using the formula from exercise 3
F = N2' * A * N1;
F = F./F(3,3);
test_detF = det(F)

%% Computes the epipolar lines l = F x1
l = F*x{1};
% Makes sure that the line has a unit normal
% makes the distance formula easier
l = l./sqrt(repmat(l(1 ,:).^2 + l(2 ,:).^2 ,[3  1]));

%% 20 points and epipolar lines in the second image
%Pick 20 points in the second image at random 
rndperm = randperm(size(x{1},2));
% plot these points in the same figure as the image
figure
imagesc(image2)
hold on
for i=1:20
    plot(x{2}(1,rndperm(i)), x{2}(2,rndperm(i)), 'ro', 'linewidth',1);
    % plot the corresponding epipolar lines using function rital.m
    rital(l(:,rndperm(i)));
end
title('CE1: 20 points with the corresponding epipolar lines')
% Are they close to each other?

%% Computes all the distances between the points
% and their corresponding epipolar lines, and plots in a histogram
figure
hist(abs(sum(l.*x{2})) ,100);
title('CE1: Distances between the points and their epipolar lines')

% The mean distance
meandist = mean(abs(sum(l.*x{2})));
%% Exercise 4
F4 = [0 1 1; 1 0 0; 0 1 1];
X1 = [1;2;3];
X2 = [3;2;1];

% epipolar e2
%e24 = svd(F4'); ??????????????
e24 = null(F4');

% Constructs the cross product matrix
e24cross = [0 -e24(3) e24(2); 
            e24(3) 0 -e24(1); 
            -e24(2) e24(1) 0];
A4 = e24cross * F4;
P24 = [A4 e24];

%% ============= Computer Exercise 2 =====================================
%  =======================================================================
%{
 Use the fundamental matrix F that you obtained in Computer Exercise 1 
to compute the camera matrices in Exercise 4. Also use triangulation 
(with DLT) to compute the 3D-points. Plot both the image, the image 
points, and the projected 3D points in the same figure. 

Don't forget to normalize when triangulating. Note that since the point 
sets are the same as in Computer Exercise 1 the normalization matrices 
will also be the same. (Alternatively one could compute cameras and 
3D points from the fundamental matrix ?F obtained with the normalized 
points and transform the cameras afterwards. This also gives a valid 
solution, but it is a different one.)
Plot the 3D-points in a 3D plot. Does it look like you expected? 
(Recall Computer Exercise 1 in Assignment 2...)
%}

% compute the camera matrices in Exercise 4
P14 = [eye(3) zeros(3,1)];

e24F = null(F');
% Constructs the cross product matrix
e24Fcross = [0 -e24F(3) e24F(2); 
            e24F(3) 0 -e24F(1); 
            -e24F(2) e24F(1) 0];
A4F = e24Fcross * F;
P24F = [A4F e24F];

%% use triangulation with DLT to compute the 3D-points
% normalizing using the normalization matrices from Computer Exercise 1
P14n = N1*P14;
P24Fn = N2*P24F;

%% Set up the DLT equations for triangulation
X = [];
for j=1:size(x{1},2)
    M2 = [P14n -[xn{1}(:,j)] [0 0 0]';
         P24Fn [0 0 0]' -[xn{2}(:,j)]];
     
    % and solve the homogeneous least squares system 
    % do this in a loop, once for each point
    [U4,S4,V4] = svd(M2);
    vstar4 = V4(:,end);
    X = [X vstar4(1:4,:)];
end

[m,n] = size(X);
homX = zeros(m,n);
for i=1:n
    homX(:,i) = pflat(X(:,i));
end

%% Project the computed points into the two images
P14X = P14*homX;
P24FX = P24F*homX;

[m4,n4] = size(P14X);
x1proj = zeros(m4,n4);
x2proj = zeros(m4,n4);
for j = 1 : n4
    x1proj(:,j) = pflat(P14X(:,j)); 
    x2proj(:,j) = pflat(P24FX(:,j));
end

%% Plot both the image, the image points, and the projected 3D points 
figure
imagesc(image1)
hold on
plot(x{1}(1,:),x{1}(2,:),'g+')
hold on;
plot(x1proj(1,:),x1proj(2,:),'ro');
title('CE2: Projected 3D points and image points into P1')
legend('Image points', 'Projected 3D points')

%% Plot the 3D-points in a 3D plot
%Does it look like you expected? 
%(Recall Computer Exercise 1 in Assignment 2...)
figure
plot3(homX(1,:), homX(2,:), homX(3,:), '.','Markersize',2)
title('CE2: The reconstructed 3D-points')
xlabel('x'), ylabel('y'), zlabel('z')


%% Useful  matlab  commands:

% Computes a 3x3 matrix containing all multiplications
% of coordinates from x1n(:,i) and x2n(:,i).
%xx = x2n(:,i)*x1n(:,i)';

% Reshapes the matrix above and adds to the M matrix
%M(i,:) = xx(:)';

% Forms an F-matrix from the solution v of the least squares problem
%Fn = reshape(v,[3  3]);

% Computes and plots all the epipolar constraints (should be roughly 0)
%plot(diag(x2n'*Fn*x1n));

% Computes the epipolar lines
%l = F*x{1};
% Makes sure that the line has a unit normal
% (makes the distance formula easier)
%l = l./sqrt(repmat(l(1 ,:).^2 + l(2 ,:).^2 ,[3  1]));

% Computes all the the distances between the points
% and  there  corresponding  lines , and  plots in a histogram
%hist(abs(sum(l.*x{2})) ,100);

% Computes the epipole
% e2 = null(F'); 
% Constructs the cross product matrix
%e2x = [0 -e2(3) e2(2); e2(3) 0 -e2(1); -e2(2) e2(1) 0]; 


