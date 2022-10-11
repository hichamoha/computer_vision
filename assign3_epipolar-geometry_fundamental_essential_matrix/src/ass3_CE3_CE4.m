
% ------------ Computer Vision - Assignment 3 -------------------------
% ------------ Hicham Mohamad, hsmo@kth.se
%{
============   3 - The Essential Matrix  ==============================

Computer Exercise 3. 
The file compEx3data.mat contains the calibration matrix K for the two 
images in Computer Exercise 1. 

Normalize the image points using the inverse of K. 

Set up the matrix M in the eight point algorithm, and solve the 
homogeneous least squares system using SVD. Check that the 
minimum singular value and Mv are both small. 

Construct the Essential matrix from the solution v. Don?t forget to make 
sure that E has two equal singular values and the third one zero. Check 
that the epipolar constraints ? x2^ E? x1 = 0 are roughly ful?lled.
%}

%% load the data
clear all, close all

load('compEx1data');
load('compEx3data');
image1 = imread('kronan1.jpg');
image2 = imread('kronan2.jpg');

%% Normalize the image points using the inverse of K
xn{1,1} = K\x{1};
xn{2,1} = K\x{2};

%% Set up the matrix M in the eight point algorithm 
% (use all the points), 
% Computes a 3x3 matrix containing all multiplications
% of coordinates from x1n(:,i) and x2n(:,i)
for i=1:size(x{1},2)
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

%% Creates a valid essential matrix from an approximate solution
Eapprox = reshape(vstar,[3 3]);

[EU,ES,EV] = svd(Eapprox);
if det(EU*EV')>0 
    E = EU*diag([1 1 0])*EV'; 
else
    EV = -EV;
    E = EU*diag([1 1 0])*EV'; 
end
test_detE = det(E)

%% Compute and plot all the epipolar constraints (should be roughly 0)
disp('The Epipolar constraints ')
xn{2}(:,1)'*E*xn{1}(:,1)
figure
plot(diag(xn{2}' * E * xn{1}),'.');
title('CE3: Plot of the epipolar constraints with E')

%% Compute the fundamental matrix
%{
Compute the fundamental matrix for the un-normalized coordinate system 
from the essential matrix and compute the epipolar lines l = Fx1. 
Pick 20 of the detected points in the second image at random and plot 
these in the same ?gure as the image. Also plot the corresponding 
epipolar lines in the same ?gure using thefunction rital.m.
%}
% un-normalized fundamental matrix from E
F = K' \ E /K;
% The essential matrix
E = E./E(3,3)

%% compute the epipolar lines l = F x1
l = F * x{1};
% Makes sure that the line has a unit normal
% makes the distance formula easier
l = l./sqrt(repmat(l(1 ,:).^2 + l(2 ,:).^2 ,[3  1]));

%% 20 points and epipolar lines in the second image
%Pick 20 points in the second image at random
%rnd = randi(length(x{1}), 20,1);
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

%% Compute the distance
%{
Compute the distance between the points and their corresponding epipolar 
lines and plot these in a histogram with 100 bins. How does this result 
compare to the corresponding result Computer Exercise 1?
%}
figure
hist(abs(sum(l.*x{2})) ,100);
title('CE3: Distances between the points and their epipolar lines')
xlabel('The different distances')
ylabel('Image points')
% The mean distance
meandist = mean(abs(sum(l.*x{2})));

%% Exercise 6
U6 = [1/sqrt(2) -1/sqrt(2) 0;
     1/sqrt(2) 1/sqrt(2) 0;
     0 0 1];
V6 = [1 0 0; 0 0 -1; 0 1 0];
E6 = U6*diag([1 1 0])*V6';

% plausible correspondence
x16 = [0; 0; 1];
x26 = [1; 1; 1];
testzero = x26'*E6*x16;

% compute s sucht that X(s) projects to x2
W6 = [0 -1 0; 1 0 0; 0 0 1];
P21 = [U6*W6*V6' U6(:,end)];
P22 = [U6*W6*V6' -U6(:,end)];
P23 = [U6*W6'*V6' U6(:,end)];
P24 = [U6*W6'*V6' -U6(:,end)];

%% =========== Computer Exercise 4 ========================================
%  ========================================================================
%{
 For the essential matrix obtained in Computer Exercise 3 compute four 
camera solutions in (10). Triangulate the points using DLT for each of 
the four camera solutions, and determine for which of the solutions the 
points are in front fo the cameras. (Since there is noise involved it 
might not be possible to find a solution with all points in front of the 
cameras. In that case select the one with the highest number of points 
in front of the cameras). 

Compute the corresponding camera matrices for the original (un-normalized) 
coordinate system and plot the image the points and the projected 
3D-points in the same figure. Does the errors look small? 

Plot the 3D points and camera centers and principal axes in a 3D plot. 
Does it look like you expected it to? 
(Recall Exercise 2 in Assignment 2...)
%}

W = [0 -1 0; 1 0 0; 0 0 1];
% Camera P1
P1 = [eye(3,3) zeros(3,1)];

% compute the four camera solutions P2
P2{1,1} = [EU*W*EV' EU(:,end)];
P2{2,1} = [EU*W*EV' -EU(:,end)];
P2{3,1} = [EU*W'*EV' EU(:,end)];
P2{4,1} = [EU*W'*EV' -EU(:,end)];

%% Triangulate the points using DLT for each of the 4 camera solutions
% Set up the DLT equations for triangulation
X = {};
for i=1:4
    Xpoints = [];
    for j=1:size(x{1},2)
        Mtriang = [P1 -[xn{1}(:,j)] [0 0 0]';
                   P2{i,1} [0 0 0]' -[xn{2}(:,j)]];
     
        % and solve the homogeneous least squares system 
        % do this in a loop, once for each point
        [EU,ES,EV] = svd(Mtriang);
        vstar = EV(:,end);
        Xpoints = [Xpoints vstar(1:4,:)];
    end
    X{i,1} = Xpoints;
end
%% homogeneous 3D-points X
[m,n] = size(X{1});
homX = {};
isInFrontOf = {};

for j=1:4
    homXpoints = zeros(m,n);
    inFrontInd = [];
    
    for i=1:n
        homXpoints(:,i) = pflat(X{j}(:,i));
        
        % determine for which of the solutions the  
        % points are in front fo the cameras
        
        % If both P1(3,:)*homX AND P2{i}(3,:)*homX are positive 
        % then the point is in front of both cameras.
        % Check if P1(3,:)*homX{j} AND P2{i}(3,:)*homX{j}.
        if (P1(3,:)*homXpoints(:,i)>0) && (P2{j}(3,:)*homXpoints(:,i)>0)
            inFrontInd = [inFrontInd i];
        end
    end
    homX{j,1} = homXpoints;
    isInFrontOf{j,1} = inFrontInd;
       
end
for i=1:4
    sz(i) = size(isInFrontOf{i},2);
end
    [~,highest] = max(sz);
disp(['The solution ', num2str(highest),+... 
      ' is with the highest number of points ' +...
      'in front of the cameras.'])

%% Project the computed points into the two images
for i=1:4
    P1n = K*P1;
    P2n{i} = K*P2{i}
    P1Xn{i,1} = P1n*homX{i};
    P2Xn{i,1} = P2n{i}*homX{i};

    [m4,n4] = size(P1Xn{1});
    proj1 = zeros(m4,n4);
    proj2 = zeros(m4,n4);
    for j = 1 : n4
        proj1(:,j) = pflat(P1Xn{i}(:,j)); 
        proj2(:,j) = pflat(P2Xn{i}(:,j));
    end
    x1proj{i,1} = proj1;
    x2proj{i,1} = proj2;
end

%% Plot both the image, the image points, and the projected 3D points 
for i=1:4
    figure
    imagesc(image2)
    hold on
    plot(x{2}(1,:),x{2}(2,:),'g+')
    hold on;
    plot(x2proj{i}(1,:),x2proj{i}(2,:),'ro');
    title(['CE4: Projected 3D points and image points into P2, solution ',num2str(i)])
    legend('Image points', 'Projected 3D points')
end

%% Plot the 3D-points in a 3D plot
%Does it look like you expected? 
%(Recall Exercise 2 in Assignment 2...)
for i=1:4
    % computer camera centers
    C{i,1} = null(P2{i});
    % compute principal axes (row 3 R3)
    viewing{i,1} = P2{i}(end, 1:3);
        
    figure
    plot3(homX{i}(1,:), homX{i}(2,:), homX{i}(3,:), '.','Markersize',2)
    hold on
    plotcams({P1; P2{i}})
    title(['CE4: The reconstructed 3D-points, solution ', num2str(i)])
    xlabel('x'), ylabel('y'), zlabel('z')
end
%% Useful matlab commands:
% Creates a valid essential matrix from an approximate solution. 
%[U,S,V] = svd(Eapprox);
%if det(U*V')>0 
%    E = U*diag([1 1 0])*V'; 
%else
%    V = -V;
%    E = U*diag([1 1 0])*V; 
%end
% Note: Computing svd on E may still give U and V that does not fulfill 
% det(U*V') = 1 since the svd is not unique. 
% So don't recompute the svd after this step.
