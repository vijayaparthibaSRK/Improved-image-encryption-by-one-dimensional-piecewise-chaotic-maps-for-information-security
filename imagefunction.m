function [imagefunction] = imagefunction( file1, file2 );
% clear all; % Erase all existing variables.


workspace; % Make sure the workspace panel is showing.
fontSize = 20;


%%%%%%%%%%%%%
% %% team: edit your file here
% %%% set the path
% folder = 'I:\F2015_486A\image';

% %reading images as array to variable 'a' & 'b'.
image1 = imread( file1 );
 image2= imread( file2 );
%
% %reading images as array to variable 'a' & 'b'.
% image1 = imread('Loaded-Mid Channel.PNG');
% image2= imread('Gross Misloaded-Mid Channel.PNG');
%


% %% info
  imfinfo( file1 )
  imfinfo( file2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% do not edit anything below here!!!!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a = image1;
b = image2;

 %% size
 a1= size(a);
 b1= size (b);

    figure;
[rows columns numberOfColorBands] = size(a);
subplot(2, 2, 1); imshow(a, []);
title('Reference image', 'Fontsize', fontSize);
set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.
redPlane = a(:, :, 1);
greenPlane = a(:, :, 2);
bluePlane = a(:, :, 3);


[rows columns numberOfColorBands] = size(b);
subplot(2, 2, 3); imshow(b, []);
title('Image2', 'Fontsize', fontSize);
set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.
redPlane1 = b(:, :, 1);
greenPlane1 = b(:, :, 2);
bluePlane1 = b(:, :, 3);

% Let's get its histograms.
[pixelCountR grayLevelsR] = imhist(redPlane);
subplot(2, 2, 2); %plot(pixelCountR, 'r');
stem(grayLevelsR,pixelCountR, 'r');
title('Histogram of Reference image', 'Fontsize', fontSize);
xlabel('Intensity values','Fontsize', fontSize);
ylabel('Number of pixels', 'Fontsize', fontSize);
xlim([0 grayLevelsR(end)]); % Scale x axis manually.

[pixelCountR1 grayLevelsR1] = imhist(redPlane1);
subplot(2, 2, 4);
%plot(pixelCountR1, 'r');
stem(grayLevelsR1,pixelCountR1);
title('Histogram of Image2', 'Fontsize', fontSize);
xlabel('Intensity values','Fontsize', fontSize);
ylabel('Number of pixels','Fontsize', fontSize);
xlim([0 grayLevelsR1(end)]); % Scale x axis manually.


%check image for different
 different = a1- b1

 if different==0
     disp('The images are same')%output display 
 else
    disp('the images are not same')
    msgbox('Test error. Please check your setting', 'Error')
    end;


