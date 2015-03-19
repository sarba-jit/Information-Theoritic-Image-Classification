QUANTIZATION_FACTOR = 32;
CELL_SIZE = 4; %greater than 4
 
for z = 1:400
        
original_img = imread( strcat('C:\Users\sarbajit\Desktop\JPEG\test\original_face\',sprintf('face%d.tif',z)),'tif' );
 
%%% MAPPER RGB -> YCbCr
%ycbcr_img = rgb2ycbcr(original_img);

ycbcr_img = original_img;
y_image = ycbcr_img(:, :, 1);
%cb_image = ycbcr_img(:, :, 2);
%cr_image = ycbcr_img(:, :, 3);
 
% rgbImage = cat(3, redImage, greenImage, blueImage);
repeat_height = size(y_image, 1)/CELL_SIZE;
repeat_width = size(y_image, 2)/CELL_SIZE;
repeat_height_mat = repmat(CELL_SIZE, [1 repeat_height]);
repeat_width_mat = repmat(CELL_SIZE, [1 repeat_width]);
y_sub_image = mat2cell(y_image, repeat_width_mat, repeat_height_mat);
%cb_sub_image = mat2cell(cb_image, repeat_width_mat, repeat_height_mat);
%cr_sub_image = mat2cell(cr_image, repeat_width_mat, repeat_height_mat);

%fprintf('Value is : %f \n',y_sub_image{8,8});
%fprintf('Value after noise : %f \n',y_sub_image{8,8}+noise(1,1,0,9.75));

%DCT TRANSFORM
for i=1:repeat_height
for j=1:repeat_width
y_sub_image{i, j} = dct2(y_sub_image{i, j});
%cb_sub_image{i, j} = dct2(cb_sub_image{i, j});
%cr_sub_image{i, j} = dct2(cr_sub_image{i, j});
end
end

%fprintf('Value before noise : %f \n',y_sub_image{8,8});
%fprintf('Value is : %f \n',y_sub_image{8,8}+noise(1,1,0,1.75));
%%Quanitzation
for i=1:repeat_height
for j=1:repeat_width
y_sub_image{i, j} = y_sub_image{i, j} / QUANTIZATION_FACTOR;
%cb_sub_image{i, j} = cb_sub_image{i, j} / QUANTIZATION_FACTOR;
%cr_sub_image{i, j} = cr_sub_image{i, j} / QUANTIZATION_FACTOR;
end
end
 
%%Entropy Encoding
for i=1:repeat_height
for j=1:repeat_width
y_sub_image{i, j}(:, 5:8) = 0;
y_sub_image{i, j}(5:8, :) = 0;
%cb_sub_image{i, j}(:, 5:8) = 0;
%cb_sub_image{i, j}(5:8, :) = 0;
%cr_sub_image{i, j}(:, 5:8) = 0;
%cr_sub_image{i, j}(5:8, :) = 0;
end
end
 
y_compressed_img = uint8(cell2mat(y_sub_image )) * QUANTIZATION_FACTOR;

%y_compressed_img = imresize(y_compressed_img, [128 128], 'bicubic');

%cb_compressed_img = uint8(cell2mat(cb_sub_image)) * QUANTIZATION_FACTOR;
%cr_compressed_img = uint8(cell2mat(cr_sub_image)) * QUANTIZATION_FACTOR;

%compressed_img = original_img;
%compressed_img(:, :, 1) = y_compressed_img;

%compressed_img(:, :, 2) = cb_compressed_img;
%compressed_img(:, :, 3) = cr_compressed_img;
%final_img_compressed = ycbcr2rgb(compressed_img);

%Noise Addition
%for i=1:128
 %   for j=1:128
  %      if (23<compressed_img(i,j)<65)
   %         r = compressed_img(i,j)*1.50;
    %        compressed_img(i,j)=compressed_img(i,j)+noise(1,1,0,r);            
     %   end 
   % end
%end


%final_img_compressed = compressed_img;

final_img_compressed = y_compressed_img;
%figure(6);
%imshow(final_img_compressed);
imwrite(final_img_compressed,strcat('C:\Users\sarbajit\Desktop\JPEG\test\4_32\random_noise\4_32_compressed\',sprintf('compressed%d.tif',z)),'tif');

%%Reconstruction of Image
 
%%Quantasize Back
%%Quanitzation
for i=1:repeat_height
for j=1:repeat_width
y_sub_image{i, j} = y_sub_image{i, j} * QUANTIZATION_FACTOR;
%cb_sub_image{i, j} = cb_sub_image{i, j} * QUANTIZATION_FACTOR;
%cr_sub_image{i, j} = cr_sub_image{i, j} * QUANTIZATION_FACTOR;
end
end
 
%%inverse dct
%DCT TRANSFORM
for i=1:repeat_height
for j=1:repeat_width
y_sub_image{i, j} = idct2(y_sub_image{i, j});
%cb_sub_image{i, j} = idct2(cb_sub_image{i, j});
%cr_sub_image{i, j} = idct2(cr_sub_image{i, j});
end
end
 
%%Stich sub images into single image
y_image = cell2mat(y_sub_image);
%cb_image = cell2mat(cb_sub_image);
%cr_image = cell2mat(cr_sub_image);
 
y_image = imresize(y_image, [128 128], 'bicubic');

%%Convert to RGB space
ycbcr_img(:, :, 1) = y_image;
%ycbcr_img(:, :, 2) = cb_image;
%ycbcr_img(:, :, 3) = cr_image;
 
%final_img = ycbcr2rgb(ycbcr_img);
 
final_img = ycbcr_img;
%figure(1);
%imshow(original_img);
 
%{
figure(2);
subplot(3, 1, 1);
title('red');
imhist(original_img(:, :, 1));
subplot(3, 1, 2);
title('green');
imhist(original_img(:, :, 2));
subplot(3, 1, 3);
title('blue');
imhist(original_img(:, :, 3));
%} 
%figure(3);
%imshow(final_img);
imwrite(final_img,strcat('C:\Users\sarbajit\Desktop\JPEG\test\4_32\random_noise\4_32_decompressed\',sprintf('decompressed%d.tif',z)),'tif');
%{ 
figure(4);
subplot(3, 1, 1);
title('Y - luma');
[counts, x] = imhist(compressed_img(:, :, 1));
stem(x, counts);
subplot(3, 1, 2);
title('Cb - Chroma Blue');
[counts, x] = imhist(compressed_img(:, :, 2));
stem(x, counts);
subplot(3, 1, 3);
title('Cr - Chroma Red');
[counts, x] = imhist(compressed_img(:, :, 3));
stem(x, counts);
%}
end