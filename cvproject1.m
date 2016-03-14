clc
clear;
myFolder = '/Users/quan199555/Documents/MATLAB/Office';
%myFolder = '/Users/quan199555/Documents/MATLAB/EnterExitCrossingPaths2cor';
%myFolder = '/Users/quan199555/Documents/MATLAB/Office';
diff1Dx = 0.5*[-1 0 1];

sigma = 4; %change sigma to get better result 
size = 6; %size of filter is n - 1 = 5
x = linspace(-size / 2, size / 2, size); 
gaussFilter = exp(-x .^ 2 / (2 * sigma ^ 2)); 
gaussFilter = gaussFilter / sum (gaussFilter); % normalize
gaussdiff1D = diff(gaussFilter);
if ~isdir(myFolder)
  errorMessage = sprintf('Error: The following folder does not exist:\n%s', myFolder);
  uiwait(warndlg(errorMessage));
  return;
end
filePattern = fullfile(myFolder, '*.jpg');
jpegFiles = dir(filePattern);
for k = 1:length(jpegFiles)
  baseFileName = jpegFiles(k).name;
  fullFileName = fullfile(myFolder, baseFileName);
  fprintf(1, 'Now reading %s\n', fullFileName);
  image = imread(fullFileName);
  greyImage = rgb2gray(image);
  sigma_n = estimate_noise(greyImage); %estimate noise
  
  greyImage = imboxfilt(greyImage, 3); %3x3 box filter
  %greyImage = imboxfilt(greyImage, 5); %5x5 box filter
  %greyImage = medfilt2(greyImage,[15 15]); %median filter
  
  %greyImage = imgaussfilt(greyImage, sigma_n); %gaussian filter
  
  ImageArray(:,:,k) = greyImage; %filtered image
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%1D diff (Does not work very well)
   % if (k > 2)
   %     imagediff = double(ImageArray(:,:,k-2)).*diff1Dx(1) + double(ImageArray(:,:,k-1)).*diff1Dx(2) + diff1Dx(3).*double(ImageArray(:,:,k));
   % else
   %     imagediff = greyImage - greyImage;
  %  end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%1D gaussian deriv (Work best + Best result)
    if (k >4)
         imagediff = double(ImageArray(:,:,k-4)).*gaussdiff1D(1) + double(ImageArray(:,:,k-3)).*gaussdiff1D(2) + gaussdiff1D(3).*double(ImageArray(:,:,k-2)) + gaussdiff1D(4).*double(ImageArray(:,:,k-1))+ gaussdiff1D(5).*double(ImageArray(:,:,k));
    else
        imagediff = greyImage - greyImage;
    end
      
  imagediff= uint8(imagediff); %convert to uint8

  level1 = graythresh(imagediff); %compute threshold
  %level(1,k) = level1; %for testing purpose
   BW1 = im2bw(imagediff,level1); %conver to binary mask

  final1 = uint8(image) .* repmat(uint8(BW1), [1,1,3]); %apply binary mask

   imshow([final1 image]);  % Display image with BW mask.
   drawnow; % Force display to update immediately.
end