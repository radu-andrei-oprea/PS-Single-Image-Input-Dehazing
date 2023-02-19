close all;clear;clc;

input = 'inputs\tree.png';
I = imread(input);

figure;
imshow(I);
title("input image");

%valori default
r = 15;
beta = 1.0; %scattering coeficient

th0 = 0.121779; th1 = 0.959710; th2 = -0.780245;
sigma = 0.041337;

tic;


%depth map
depth = getDepth(I, th0, th1, th2, sigma);

%lumina atmosferica
A = estimateAtmLight(I, depth);

%transmisie
transm = transmission(beta, depth);

J = radiance(I, A, transm);

toc;

%figures
figure;
imshow(depth);
title("depth map");

figure;
imshow(transm);
title("transmission");

figure;
imshow(J);
title("dehazed image");

function [depthMap] = getDepth(I, th0, th1, th2, sigma)
   
   %hue, saturation, value map
   hsvmap = rgb2hsv(I);
    
   %saturation        value
   s = hsvmap(:,:,2); v = hsvmap(:,:,3);
   
   %random array
   mat = normrnd(0, sigma, size(I, 1), size(I, 2));
    
   %formula
   depthMap = th0 + th1*v + th2*s + mat;

   imwrite(depthMap, "outputs\depthMap.jpg");
end

function [A] = estimateAtmLight(I, dm)
    
    height = size(I, 1); width = size(I, 2);
    I = double(I)/255;
    
    %nr de 0.1% cei mai luminosi pixeli
    noOfBrightest= ceil(0.001 * height * width);
    
    %depthmap sortat, locatia originala a pixelilor
    [sorted, location] = sort(dm(:));

    %rearanjarea imaginii originale in coloane
    realign = reshape(I, height * width, 1, 3);
    
    %array pentru pixelii cu intensitate cea mai mare
    int_arr = zeros(noOfBrightest, 1, 3);

    %array pentru normele pixelilor
    norm_arr = zeros(noOfBrightest, 1);

    %popularea celor doua array-uri
    for i = 1: noOfBrightest
        p = location(height * width + 1 - i);
        int_arr(i, 1, :) = realign(p, 1, :);
        norm_arr(i) = norm(int_arr(i,:));
    end

    [sorted, location] = sort(norm_arr(:));
    A = int_arr(location(noOfBrightest - length(location) + 1 : noOfBrightest),:);
    
    A = max(A);
end

function [transm] = transmission(beta, depth)
    
    %formula
    transm = exp(-beta*depth);
    imwrite(transm, "outputs\transmission.jpg");
end

function [J] = radiance(I, A, transm)
    
    I = double(I)/255;
    
    %formula: J = (I - A) / tr + A
    % tr : [0.1, 0.9]

    height = size(transm, 1); width = size(transm, 2);
   
    %reducere noise
    tmin = 0.1;
    tmax = 0.9;
    
    %radianta
    J = zeros(height, width);
     
    for i = 1:3
        J(:,:,i) = I(:,:,i)-A(i);
    end
    
    t = transm;
    %corectare transmisie
    for i = 1:height
        for j = 1:width
            if t(i,j) > tmax
                t(i,j) = tmax;
            elseif t(i,j) < tmin
                t(i,j) = tmin;
            end
        end
    end

    for i = 1:3
        J(:,:,i) = J(:,:,i)./t;
        J(:,:,i) = J(:,:,i) + A(i);
    end

    imwrite(J, "outputs\radiance.jpg");
    
end
