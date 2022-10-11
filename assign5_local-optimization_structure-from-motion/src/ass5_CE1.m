
% ------------ Computer Vision - Assignment 5 -------------------------
% ------------ Hicham Mohamad, hsmo@kth.se

% Computer Exercise 1
%{
 In Computer Exercise 3 and 4 of Assignment 3, you computed a solution to 
the two-view structure form motion problem for the two images of Figure 1 
using the 8-point algorithm. In this exercise the goal is to use the 
solution from Assignment3 as a starting solution and locally improve 
it using the Levenberg-Marquardt method.

The file LinearizeReprojErr.m contains a function that for a given set of 
cameras, 3D points and imagepoints, computes the linearization (6). 
The file update_solution.m contains a function that computes a new set of 
cameras and 3D points from an update deltav computed by any method. 
The file ComputeReprojectionError.m computes the reprojection error for a 
given set of cameras, 3D points and image points. It also returns the 
values of all the individual residuals as a second output. 

In the Levenberg-Maquardt method the update is given by 
deltav = -(J(vk)^T J(vk) + lambdaI)^{-1} J(vk)^T r(vk). (12) 

Using this scheme and starting from the solution that you got in 
Assignment 3, plot the reprojection error versus the iteration number 
for lambda = 1. Also plot histograms of all the residual values before and 
after running the Levenberg-Maquardt method. 

Try varying lambda. What happens if lambda is very large/small?

%}
%% running In Computer Exercise 3 and 4 of Assignment3
clear all, close all
ass3_CE3_CE4

%% use the solution from Assignment3 as a starting solution
close all

% starting from the solution obtained in Assignment 3
% set of cameras
P = {P1,P2{highest}};
% 3D points
U = homX{highest};
% image points
u = xn;

lambda = 1e-15;
iter = 25;

%% Compute the reprojection error 
startErr = zeros(iter,1);
for i=1:iter
    % from the current solution P,U and the imagedata u. 
    % The value of each residual is in res.
    [startErr(i),startRes] = ComputeReprojectionError(P,U,u);
end
%% plot the reprojection error versus the iteration number for lambda = 1.
figure
plot(1:iter, startErr)
%semilogy(1:iter, startErr)
title('Reprojection error before Levenberg, \lambda=1')
xlabel('Iteration number')

%% plot histograms of all the residual values 
% before and after running the Levenberg-Maquardt method
figure
histogram(startRes, 100)
title('Residual values before Lavenberg-Maquardt, \lambda=1')

%% improve the solution using the Levenberg-Marquardt method
err = zeros(iter,1);
Pnew = P;
Unew = U;
for i=1:iter
    %Computes the reprejection error and the values of all the residuals 
    %for the current solution P,U,u. 
    [err(i),res] = ComputeReprojectionError(Pnew,Unew,u);

    %Computes the r and J matrices for the appoximate 
    %linear least squares problem. 
    [r,J] = LinearizeReprojErr(Pnew,Unew,u);

    % Computes the LM update. 
    C = J'*J+lambda*speye(size(J,2)); 
    c = J'*r; 
    deltav = -C\c;

    %Updates the variabels 
    [Pnew ,Unew] = update_solution(deltav ,Pnew,Unew);
end
    

%% plot histograms of all the residual values 
% after running the Levenberg-Maquardt method
figure
histogram(res, 100)
title('Residual values after Lavenberg-Maquardt, \lambda=1e-15')

%% plot the reprojection error versus the iteration number for lambda = 1.
figure
plot(1:iter, err)
title('Reprojection error after Levenberg, \lambda=1e-15')
xlabel('the iteration number')

%% Useful matlab commands:
%{
%Takes two camera matrices and puts them in a cell. 
P = {P1,P2}

%Computes the reprejection error and the values of all the residuals 
%for the current solution P,U,u. 
[err,res] = ComputeReprojectionError(P,U,u);

%Computes the r and J matrices for the appoximate 
%linear least squares problem. 
[r,J] = LinearizeReprojErr(P,U,u)

% Computes the LM update. 
C = J'*J+lambda*speye(size(J,2)); 
c = J'*r; 
deltav = -C\c;

%Updates the variabels 
[Pnew ,Unew] = update_solution(deltav ,P,U);

%}