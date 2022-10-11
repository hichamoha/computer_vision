
% ------------ Computer Vision - Assignment 4 -------------------------
% ------------ Hicham Mohamad, hsmo@kth.se    -------------------------
%{
============   2 - Plane Fitting  ==================================

Computer Exercise 1. 
Consider two images of a house and a set of 3D points from 
the walls of the house. The goal of this exercise is to estimate 
the location of the wall with the most 3D points. The ?le 
compEx1data.mat contains cameras P, inner parameters K for both 
cameras, scene points X and some extra points x from image 1. 
%}

%% load the data
clear all, close all
load('compEx1data');
image1 = imread('house1.jpg');
image2 = imread('house2.jpg');

%% Solve the total least squares problem with all the points
%{
Solve the total least squares problem with all the points. 
Compute the RMS distance between the 3D-points and the plane.
%}
[m,n] = size(X);
homX = zeros(m,n);
for i=1:n
    homX(:,i) = pflat(X(:,i));
end

%Computes the mean 3D point 
meanX = mean(homX,2); 
%Subtracts the mean from the 3D points 
Xtilde = (homX - repmat(meanX ,[1 size(homX,2)])); 
%Computes the matrix from Exercise 2 
M = Xtilde(1:3,:)*Xtilde(1:3,:)'; 
%Computes eigenvalues and eigenvectors of M
[V,D] = eig(M); 

% Extract the eigenvector that correspond to a,b,c
% to the smallest eigenvalue which is in the first column
abc = V(:,1);
% compute the optimal d after taking the derivative
d = - abc'*meanX(1:3);

%Computes a plane from the obtained parameters a,b,c,d.
planeParams = [abc; d]; 

%Makes sure that the plane has a unit length norm
plane = planeParams./norm(planeParams(1:3));

%Computes the RMS error
RMS = sqrt(sum((plane'*homX).^2)/size(homX,2));
disp(['RMS distance ', num2str(RMS)])

%% Fitting the plane using RANSAC
%{
Use RANSAC to robustly ?t a plane to the 3D points X. If a 3D point 
is an inlier when its distance to the plane is less than 0.1, how many 
inliers do you get? Compute the RMS distance between the plane 
obtained with RANSAC and the distance to the 3D points. 
Is there any improvement? Plot the absolute distances between the 
plane and the points in a histogram with 100 bins.
%}
iters = 10;
RANSACplane = {};
inliersNbr = zeros(iters,1);
inliersIndx = {iters; 1};
for i=1:iters
    % pseudo random integers
    randind = randi(n, [3 1]);
    %Computes a plane from a sample set.
    RANSACplane{i,1} = null(homX(:,randind)'); 

    %Makes sure that the plane has a unit length norm
    RANSACplane{i} = RANSACplane{i}./norm(RANSACplane{i}(1:3));

    indexSet = [];
    for j=1:n
        % Finds the the indices for which the distance to 
        % the plane is less than 0.1. 
        % Note: Works only if the 4th coordinate of all the 
        % points in X is 1.
        if(abs(RANSACplane{i}'*homX(:,j)) < 0.1)
            inliersNbr(i) = inliersNbr(i) + 1 ; 
            indexSet = [indexSet j];
        end
    end
    inliersIndx{i} = indexSet;
end 

%%  select the solution that gives the largest consensus set
largestConsensusIdx = find(inliersNbr == max(inliersNbr));
consensusSize = inliersNbr(largestConsensusIdx);
disp(['Number of inliers with RANSAC ',num2str(consensusSize)])

bestPlane = RANSACplane{largestConsensusIdx};
bestX = homX(:, inliersIndx{largestConsensusIdx});

% Computes the RMS error
ransacRMS = sqrt(sum((bestPlane'*bestX).^2)/size(bestX,2));
disp(['RMS distance with RANSAC ', num2str(ransacRMS)])
 
figure
hist(abs(bestPlane'*bestX) ,100);
title('CE1: Absolute distances using RANSAC')
xlabel('Different distances')
ylabel('All the 3D points')

%% Solving the total least squares problem with only the inliers
%{
Solve the total least squares problem with only the inliers. 
Compute the RMS distance between the 3D-points and the new plane. 
Plot the absolute distances between the plane and the points in 
a histogram with 100 bins. Which estimation produced the best result? 
%}
%Computes the mean 3D point 
bestmeanX = mean(bestX,2); 
%Subtracts the mean from the 3D points 
bestXtilde = (bestX - repmat(bestmeanX ,[1 size(bestX,2)])); 
%Computes the matrix from Exercise 2 
bestM = bestXtilde(1:3,:)*bestXtilde(1:3,:)'; 
%Computes eigenvalues and eigenvectors of M
[bestV,bestD] = eig(bestM); 

% Extract the eigenvector that correspond to a,b,c
% to the smallest eigenvalue which is in the first column
bestabc = bestV(:,1);
% compute the optimal d after taking the derivative
bestd = - bestabc'*bestmeanX(1:3);

%Computes a plane from the obtained parameters a,b,c,d.
bestplaneParams = [bestabc; bestd]; 

%Makes sure that the plane has a unit length norm
inliersPlane = bestplaneParams./norm(bestplaneParams(1:3));

%Computes the RMS error
inliersRMS = sqrt(sum((inliersPlane'*bestX).^2)/size(bestX,2));
disp(['Inliers RMS distance ', num2str(inliersRMS)])

figure
hist(abs(inliersPlane'*bestX) ,100);
title('CE1: TLS absolute distances with only inliers')
xlabel('Different distances')
ylabel('Inlier 3D points')

%% Plot the projection of the inliers into the images. 
% Where are these located? 
inliersProj = {2;1};

for i=1:2
    PbestX{i,1} = P{i}*bestX;
    for j=1:size(bestX,2)
        bestXproj(:,j) = pflat(PbestX{i}(:,j));
    end
    inliersProj{i,1} = bestXproj;
end

figure
imagesc(image1)
hold on
plot(inliersProj{1}(1,:),inliersProj{1}(2,:),'r.');
title('CE1: Projected inliers into image 1')
%legend('Image points', 'Projected 3D points')

figure
imagesc(image2)
hold on
plot(inliersProj{2}(1,:),inliersProj{2}(2,:),'r.');
title('CE1: Projected inliers into image 2')
%legend('Image points', 'Projected 3D points')
    
%% Homography from camera 1 to camera 2
%{
Using the method in assignment 1, (Exercise 6,) compute a 
homography from camera 1 to camera 2. (Don?t forget that the formula 
only works for normalized cameras.) Plot the points x in image 1. 
Transform the points using the homography and plot them in image 2. 
Which ones seem to be correct, and why?
%}
% normalize the cameras
Pn{1} = K\P{1};
Pn{2} = K\P{2};

% Rotation and translation
R = Pn{2}(1:3,1:3);
t = Pn{2}(:,4);

homoPi = pflat(inliersPlane);
pi = homoPi(1:3);
% compute the homography from camera1 to camera2: 
H = (R - t*pi');

%% Plot the points x in image 1
figure
imagesc(image1)
hold on
plot(x(1,:),x(2,:),'r*');
title('CE1: The extra points x in image 1')

%% Transform the points using the homography and plot them in image 2
x2 = H*x;
figure
imagesc(image2)
hold on
plot(x2(1,:),x2(2,:),'g*');
title('CE1: The transformed points x in image 2')

%% Useful matlab commands:
%{
%Computes the mean 3D point 
meanX = mean(X,2); 
%Subtracts the mean from the 3D points 
Xtilde = (X - repmat(meanX ,[1 size(X,2)])); 
%Computes the matrix from Exercise 2 
M = Xtilde(1:3,:)*Xtilde(1:3,:)'; 
%Computes eigenvalues and eigenvectors of M
[V,D] = eig(M); 

%Computes a plane from a sample set.
plane = null(X(:,randind)'); 

%Makes sure that the plane has a unit length norm
plane = plane./norm(plane(1:3));

% Finds the the indices for which the distance to 
% the plane is less than 0.1. 
%Note: Works only if the 4th coordinate of all the points in X is 1.
inliers = abs(plane'*X) <= 0.1; 

%Computes the RMS error
RMS = sqrt(sum((plane'*X).^2)/size(X,2)); 
%}
