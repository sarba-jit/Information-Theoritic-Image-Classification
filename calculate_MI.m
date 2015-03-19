%class 1
I1 = imread('C:\Users\sarbajit\Desktop\JPEG\test\4_32\random_noise\test\face1.tif');
I2 = imread('C:\Users\sarbajit\Desktop\JPEG\test\4_32\random_noise\test\compressed1.tif');
I3 = imread('C:\Users\sarbajit\Desktop\JPEG\test\4_32\random_noise\test\decompressed1.tif');
I4 = imread('C:\Users\sarbajit\Desktop\JPEG\test\4_32\random_noise\test\face2.tif' );
I5 = imread('C:\Users\sarbajit\Desktop\JPEG\test\4_32\random_noise\test\compressed2.tif' );
I6 = imread('C:\Users\sarbajit\Desktop\JPEG\test\4_32\random_noise\test\decompressed2.tif' );
I7 = imread('C:\Users\sarbajit\Desktop\JPEG\test\4_32\random_noise\test\face3.tif' );
I8 = imread('C:\Users\sarbajit\Desktop\JPEG\test\4_32\random_noise\test\compressed3.tif' );
I9 = imread('C:\Users\sarbajit\Desktop\JPEG\test\4_32\random_noise\test\decompressed3.tif' );

%class 2
J1 = imread('C:\Users\sarbajit\Desktop\JPEG\test\4_32\random_noise\test\face11.tif');
J2 = imread('C:\Users\sarbajit\Desktop\JPEG\test\4_32\random_noise\test\compressed11.tif');
J3 = imread('C:\Users\sarbajit\Desktop\JPEG\test\4_32\random_noise\test\decompressed11.tif');
J4 = imread('C:\Users\sarbajit\Desktop\JPEG\test\4_32\random_noise\test\face12.tif' );
J5 = imread('C:\Users\sarbajit\Desktop\JPEG\test\4_32\random_noise\test\compressed12.tif' );
J6 = imread('C:\Users\sarbajit\Desktop\JPEG\test\4_32\random_noise\test\decompressed12.tif' );
J7 = imread('C:\Users\sarbajit\Desktop\JPEG\test\4_32\random_noise\test\face13.tif' );
J8 = imread('C:\Users\sarbajit\Desktop\JPEG\test\4_32\random_noise\test\compressed13.tif' );
J9 = imread('C:\Users\sarbajit\Desktop\JPEG\test\4_32\random_noise\test\decompressed13.tif' );

%class 3
K1 = imread('C:\Users\sarbajit\Desktop\JPEG\test\4_32\random_noise\test\face23.tif');
K2 = imread('C:\Users\sarbajit\Desktop\JPEG\test\4_32\random_noise\test\face24.tif' );
K3 = imread('C:\Users\sarbajit\Desktop\JPEG\test\4_32\random_noise\test\face25.tif' );
K4 = imread('C:\Users\sarbajit\Desktop\JPEG\test\4_32\random_noise\test\compressed25.tif' );
K5 = imread('C:\Users\sarbajit\Desktop\JPEG\test\4_32\random_noise\test\decompressed25.tif' );
K6 = imread('C:\Users\sarbajit\Desktop\JPEG\test\4_32\random_noise\test\face26.tif' );

%class 4
L1 = imread('C:\Users\sarbajit\Desktop\JPEG\test\4_32\random_noise\test\face37.tif');
L4 = imread('C:\Users\sarbajit\Desktop\JPEG\test\4_32\random_noise\test\compressed37.tif' );
L5 = imread('C:\Users\sarbajit\Desktop\JPEG\test\4_32\random_noise\test\decompressed37.tif' );

%class 5
M1 = imread('C:\Users\sarbajit\Desktop\JPEG\test\4_32\random_noise\test\face47.tif');
M2 = imread('C:\Users\sarbajit\Desktop\JPEG\test\4_32\random_noise\test\face43.tif' );
M3 = imread('C:\Users\sarbajit\Desktop\JPEG\test\4_32\random_noise\test\face42.tif' );
M4 = imread('C:\Users\sarbajit\Desktop\JPEG\test\4_32\random_noise\test\compressed42.tif' );
M5 = imread('C:\Users\sarbajit\Desktop\JPEG\test\4_32\random_noise\test\decompressed42.tif' );
M6 = imread('C:\Users\sarbajit\Desktop\JPEG\test\4_32\random_noise\test\face49.tif' );

A= mutualInfo(J4,J3);
B= mutualInfo(J4,J1);
D = mutualInfo(J4,J4);
%B = mutalInfo_MIT(I6,I7);
C = entropy(J4);
%D = conditionalEntropy(K6,K7);


fprintf('Mutal Info  or-dec is %.5f. \n', A);
fprintf('Mutal Info  or-com is %.5f. \n', B);
fprintf('Mutal Info  or-or %.5f. \n', D);
fprintf('Entropy is %.5f. \n',C);
