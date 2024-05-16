
clc;
clear all;
close all;
% Demo by Image Analyst.
  % Close all figures (except those of imtool.)clear;clc;

img = imread('lena.jpg');
img = rgb2gray(img);
figure(1)
imshow(img);
title('Original image ');
tic
timg = img;
r = 3.62;
x(1) = 0.7;
row = size(img,1);
col = size(img,2);
s = row*col;
%Creation of Logistic function
for n=1:s-1
    x(n+1) = r*x(n)*(1-x(n));
end

[so,in] = sort(x);

%Start of Confusion
timg = timg(:);
for m = 1:size(timg,1)
    
    t1 = timg(m);
    timg(m)=timg(in(m));
    timg(in(m))=t1;
    
end
%End of confussion


%Creation of diffusion key

p=3.628;
k(1)=0.632;
for n=1:s-1
    k(n+1) = cos(p*acos(k(n)));
end
k = abs(round(k*255));

ktemp = de2bi(k);
ktemp = circshift(ktemp,1);
ktemp = bi2de(ktemp)';
key = bitxor(k,ktemp);

%Ending creation of diffusion key


%Final Encryption Starts
timg = timg';
timg = bitxor(uint8(key),uint8(timg));
himg = reshape(timg,[row col]);
figure(9)
imshow(himg);
title('Encryption image ');
%Entropy Test
clearvars;
workspace;  % Make sure the workspace panel is showing.
format long g;
format compact;
fontSize = 16;
numberOfIterations = 1200;
K = 1.3;
p = ones(numberOfIterations, 1);
theta = zeros(numberOfIterations, 1);
plotColors = jet(numberOfIterations);
originalImage = imread('lena.jpg');
outputImage = zeros(size(originalImage), 'uint8');
[rows, columns, numberOfColorchannels] = size(originalImage)
for col = 1 : columns
    fprintf('Column %d of %d.\n', col, columns)
    for row = 1 : rows
        % Find initial p and theta for this pixel.
        p1 = sqrt((row-rows/2).^2 + (col-columns/2).^2);
        theta1 = atan((row-rows/2)/(col-columns/2));
        % Find next p and theta for this pixel.  
        % Where are we going to send this pixel to?
        % Compute new p according to the formula.
        p2 = p1 + K * sin(theta1);
        % Compute new theta according to the formula.
        theta2 = theta1 + p1;
        % Make them mod 2 * pi
        p2 = mod(p2, 2*pi);
        theta2 = mod(theta2, 2*pi);
        % Convert from polar to cartesian coordinates.
        row2 = round(p2 * sin(theta2) + rows/2);
        col2 = round(p2 * cos(theta2) + columns/2);
        if isnan(row2) || isnan(col2)
            continue;
        end
        % Skip it if it's outside the image.
        if row2 <= 1 || row2 > rows || col2 <= 1 || col2 > columns
            continue
        end
        outputImage(row2, col2, :) = originalImage(row, col, :);
    end
end
% Display output image
figure(2);
imshow(outputImage)
title('Complete Chaos!', 'FontSize', fontSize);
axis('on', 'image');
impixelinfo
img = imread('1.jpg');
img = rgb2gray(img);
timg = img;
r = 3.62;
x(1) = 0.7;
row = size(img,1);
col = size(img,2);
s = row*col;
%Creation of Logistic function
for n=1:s-1
    x(n+1) = r*x(n)*(1-x(n));
end

[so,in] = sort(x);

%Start of Confusion
timg = timg(:);
for m = 1:size(timg,1)
    
    t1 = timg(m);
    timg(m)=timg(in(m));
    timg(in(m))=t1;
    
end
%End of confussion


%Creation of diffusion key

p=3.628;
k(1)=0.632;
for n=1:s-1
    k(n+1) = cos(p*acos(k(n)));
end
k = abs(round(k*255));

ktemp = de2bi(k);
ktemp = circshift(ktemp,1);
ktemp = bi2de(ktemp)';
key = bitxor(k,ktemp);

%Ending creation of diffusion key


%Final Encryption Starts
timg = timg';
timg = bitxor(uint8(key),uint8(timg));
himg = reshape(timg,[row col]);

%Final Encryption Ends
toc
%Decryption Start
timg = bitxor(uint8(key),uint8(timg));
timg = timg(:);
for m = size(timg,1):-1:1
    
    t1 = timg(m);
    timg(m)=timg(in(m));
    timg(in(m))=t1;
    
end
%Decryption End
timg = reshape(timg,[row col]);


figure(3)
imshow(timg);
title('Decryption image ');

figure(4)
imhist(img);
title('Original image histogram');

figure(5)
imhist(himg);
title('Decrypted image histogram');
% Reading all three images
F = imread('1.jpg');
S = imread('2.jpg');
V = imread('lena.jpg');

%Converting color images to Grayscale
F = im2double(rgb2gray(F));
S = im2double(rgb2gray(S));
V = im2double(rgb2gray(V));

[rows cols] = size(F);
Tmp = [];
Tmp1 = [];
temp = 0;

% Saving the patch(rows x 5 columns) of second(S) & third(V) images in    
% S1 & V1 resp for future use.
for i = 1:rows
    for j = 1:5
        S1(i,j) = S(i,j);
        V1(i,j) = V(i,j);
    end
end

% Performing Correlation i.e. Comparing the (rows x 5 column) patch of 
% first image with patch of second image i.e. S1 saved earlier.
for k = 0:cols-5 % (cols - 5) prevents j from going beyond boundary of image.
    for j = 1:5
        F1(:,j) = F(:,k+j);% Forming patch of rows x 5 each time till cols-5.
    end
    temp = corr2(F1,S1);% comparing the patches using correlation.
    Tmp = [Tmp temp]; % Tmp keeps growing, forming a matrix of 1*cols
    temp = 0;
end

[Min_value, Index] = max(Tmp);% Gets the Index with maximum value from Tmp.

% Determining the number of columns of new image. Rows remain the same.
n_cols = Index + cols - 1;

Opimg = [];
for i = 1:rows
    for j = 1:Index-1
        Opimg(i,j) = F(i,j);% First image is pasted till Index.
    end
    for k = Index:n_cols
        Opimg(i,k) = S(i,k-Index+1);%Second image is pasted after Index.
    end    
end

[r_Opimg c_Opimg] = size(Opimg);

% Performing Correlation i.e. Comparing the (rows x 5 column) patch of 
% second image with patch of third image i.e. V1 saved earlier.
for k = 0:c_Opimg-5% to prevent j to go beyond boundaries.
    for j = 1:5
        Opimg1(:,j) = Opimg(:,k+j);% Forming patch of rows x 5 each time till cols-5.
    end
    temp = corr2(Opimg1,V1);% comparing the patches using correlation.
    Tmp1 = [Tmp1 temp]; % Tmp keeps growing, forming a matrix of 1*cols
    temp = 0;
end

% Determining the size of third image for future use.
[r_V, c_V] = size(V);

[Min_value, Index] = max(Tmp1);

% Determining new column for final stitched image.
% Rows remain the same.
n_cols = Index + c_V - 1;

Opimg1 = [];
for i = 1:rows
    for j = 1:Index-1
        Opimg1(i,j) = Opimg(i,j);% Previous stitched image is pasted till new Index.
    end
    for k = Index:n_cols
        Opimg1(i,k) = V(i,k-Index+1);%Third image is pasted after new Index.
    end    
end
% Determining the size of Final Stitched image.
[r_Opimg c_Opimg] = size(Opimg1);

figure,
subplot(2,3,1);
imshow(F);axis ([1 c_Opimg 1 r_Opimg])
title('First Image');

subplot(2,3,2);
imshow(S);axis ([1 c_Opimg 1 r_Opimg])
title('Second Image');

subplot(2,3,3);
imshow(V);axis ([1 c_Opimg 1 r_Opimg])
title('Third Image');

subplot(2,3,[4 5 6]);% Final Stitched image should get most of the space in subplot.
imshow(Opimg1);axis ([1 c_Opimg 1 r_Opimg])
title('Decrypted Image');
imagefunction('2.jpg','3.jpg')
% End of Code