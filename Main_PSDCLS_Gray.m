clc
clear
close all
format compact

%% Plain Image
PlainImg=imread('lena_gray_256.tif');
%PlainImg=imread('Tank256.bmp');
%PlainImg=imread('lena_gray_512.tif');
%PlainImg=imread('Boat512.bmp');
%PlainImg=imread('lena_gray_1024.tif');
%PlainImg=imread('Cameraman1024.bmp');

PlainImg=double(PlainImg);
N=size(PlainImg,1);

%% Generate Key
Key = '6b679b3c77826d30a79e612114a8c18df984c176f4e529f684748ad052241b17';
H1=HashFunction(PlainImg,'SHA-256');
%%  Encryption and Decryption
Runsize=1;
for i=1:Runsize
    tic
        [EncImg]=Encryption(PlainImg,Key,H1);
    tEnc(i)=toc;
    tic
        DecImg=Decryption(EncImg,Key,H1);
    tDec(i)=toc;
end

fprintf('%d Runs --> Average Encryption Time = %f\n',Runsize,sum(tEnc)/Runsize);
fprintf('%d Runs --> Average Decryption Time = %f',Runsize,sum(tDec)/Runsize);
%%
figure('Name','Simulation Result ...')
subplot(2,3,1),imshow(uint8(PlainImg)),  title('Plain Image');
subplot(2,3,2),imshow(uint8(EncImg)),    title('Cipher Image');
subplot(2,3,3),imshow(uint8(DecImg)),    title('Decrypted Image');
subplot(2,3,4),imhist(uint8(PlainImg)),  title('Histogram of Plain Image');
subplot(2,3,5),imhist(uint8(EncImg)),    title('Histogram of Cipher Image');
subplot(2,3,6),imhist(uint8(DecImg)),    title('Histogram of Decrypted Image');
% Test Parameter
DifPlainDec=sum(abs(uint8(PlainImg(:))-uint8(DecImg(:))));
fprintf('\n\n|PlainImage - DecryptedImage| = %d',DifPlainDec);

% Chi-square test results of cipher images.
[VarP, ChiP, VarC, ChiC] =Var_Chisquare_Gray(PlainImg,EncImg);
fprintf('\ncipher images variance = %d  Chi-square = %d',VarC, ChiC);
% Entropy
PlainImg_Entropy = Entropy(PlainImg,N);
fprintf('\n\nPlainImage  Entropy = %f',PlainImg_Entropy);
EncImg_Entropy = Entropy(EncImg,N);
fprintf('\nCipherImage Entropy = %f',EncImg_Entropy);

% Correlation Coefficient
fprintf('\n\nCorrelation Coefficient:\n');
CC=AdjancyCorrPixelRandNew(PlainImg,EncImg);
disp(CC);

% MSE
MSEPE=mse(uint8(EncImg),uint8(PlainImg));
R=msefun(EncImg,PlainImg);
fprintf('\n MSE PlainImage and Encrypted Image: %f  %f', MSEPE,R);
MSEPD=mse(uint8(DecImg),uint8(PlainImg));
R=msefun(DecImg,PlainImg);
fprintf('\n MSE PlainImage and Decrypted Image: %f  %f', MSEPD , R);

% PSNR
PSNRCPE=psnr(uint8(EncImg),uint8(PlainImg));
fprintf('\n PSNR PlainImage and Encrypted Image: %f ', PSNRCPE);
PSNRCPD=psnr(uint8(DecImg),uint8(PlainImg));
fprintf('\n PSNR PlainImage and Decrypted Image: %f ', PSNRCPD);

% NC
normcoPlainEnc=normcor(uint8(EncImg),uint8(PlainImg));
fprintf('\n Normalized cross-correlation function (NC) Plain-Cipher Image: %f', normcoPlainEnc);
normcoPlainDec=normcor(uint8(DecImg),uint8(PlainImg));
fprintf('\n Normalized cross-correlation function (NC) Plain-Decrypted Image: %f', normcoPlainDec);

% SSIM
[ssimvalPE] = ssim(uint8(EncImg),uint8(PlainImg));
fprintf('\n SSIM PlainImage and Encrypted Image: %f ', ssimvalPE);
[ssimvalPD] = ssim(uint8(DecImg),uint8(PlainImg));
fprintf('\n SSIM PlainImage and Decrypted Image: %f ', ssimvalPD);




%% NPCR and UACI Test
%Encryption and Decryption for 1 bit change in Plain Image
fprintf('\n\n NPCR and UACI Test ...')
for i=1:1
    PlainImg1bit=PlainImg;      
    pos1=1+floor(rand(1)*N);
    pos2=1+floor(rand(1)*N);
    
    fprintf('\nBefore change 1 bit of PlainImage at location (%d,%d) = %d',pos1,pos2,PlainImg1bit(pos1,pos2));
    PlainImg1bit(pos1,pos2) =uint8(mod(double(PlainImg1bit(pos1,pos2)+1),256));
    fprintf('\nAfter change 1 bit of PlainImage at location (%d,%d) = %d',pos1,pos2,PlainImg1bit(pos1,pos2));
    H1bit=HashFunction(PlainImg1bit,'SHA-256');
    [EncImg1bit]=Encryption(PlainImg1bit,Key,H1bit);
    [DecImg1bit] = Decryption(EncImg1bit,Key,H1bit);
    
    [npcr1, uaci1]= NPCR_UACI(uint8(EncImg),uint8(EncImg1bit));
    fprintf('\nNPCR = %f   UACI=%f \n',npcr1, uaci1);
end
%%  Cropping attack
fprintf('\nCropping attack ... \n')
figure('Name','Cropping attack ...'),
if N==256
    crop64=32;
    crop16=64;
    crop4=128;
    crop2=256;
elseif N==512
    crop64=64;
    crop16=128;
    crop4=256;
    crop2=512;
elseif N==1024
    crop64=128;
    crop16=256;
    crop4=512;
    crop2=1024;
end
%*****************************************************************************
EncImgCrop64=EncImg;
EncImgCrop64(1:crop64, 1:crop64)=0;
DecImgCrop64 = Decryption(EncImgCrop64,Key,H1);
PSNRCrop64=psnr(uint8(DecImgCrop64),uint8(PlainImg));
[ssimvalCrop64] = ssim(uint8(DecImgCrop64),uint8(PlainImg));
fprintf('1/64  cropped cipher image PSNR= %f  SSIM= %f \n',PSNRCrop64, ssimvalCrop64);
subplot(2,4,1), imshow(uint8(EncImgCrop64)),  title('1/64 Cropped Cipher Image ');
subplot(2,4,5), imshow(uint8(DecImgCrop64)),  title('Decrypted Image ');
%*****************************************************************************
EncImgCrop16=EncImg;
EncImgCrop16(1:crop16, 1:crop16)=0;
DecImgCrop16 = Decryption(EncImgCrop16,Key,H1);
PSNRCrop16=psnr(uint8(DecImgCrop16),uint8(PlainImg));
[ssimvalCrop16] = ssim(uint8(DecImgCrop16),uint8(PlainImg));
fprintf('1/16  cropped cipher image PSNR= %f  SSIM= %f \n',PSNRCrop16, ssimvalCrop16);
subplot(2,4,2), imshow(uint8(EncImgCrop16)),  title('1/16 Cropped Cipher Image ');
subplot(2,4,6), imshow(uint8(DecImgCrop16)),  title('Decrypted Image ');
%*****************************************************************************************
EncImgCrop4=EncImg;
EncImgCrop4(1:crop4, 1:crop4)=0;
DecImgCrop4 = Decryption(EncImgCrop4,Key,H1);
PSNRCrop4=psnr(uint8(DecImgCrop4),uint8(PlainImg));
[ssimvalCrop4] = ssim(uint8(DecImgCrop4),uint8(PlainImg));
fprintf('1/4   cropped cipher image PSNR= %f  SSIM= %f \n',PSNRCrop4, ssimvalCrop4);
subplot(2,4,3), imshow(uint8(EncImgCrop4)),  title('1/4 Cropped Cipher Image ');
subplot(2,4,7), imshow(uint8(DecImgCrop4)),  title('Decrypted Image ');
%*****************************************************************************************
EncImgCrop2=EncImg;
EncImgCrop2(1:crop4, 1:crop2)=0;
DecImgCrop2 = Decryption(EncImgCrop2,Key,H1);
PSNRCrop2=psnr(uint8(DecImgCrop2),uint8(PlainImg));
[ssimvalCrop2] = ssim(uint8(DecImgCrop2),uint8(PlainImg));
fprintf('1/2   cropped cipher image PSNR= %f  SSIM= %f \n',PSNRCrop2, ssimvalCrop2);
subplot(2,4,4), imshow(uint8(EncImgCrop2)),  title('1/2 Cropped Cipher Image ');
subplot(2,4,8), imshow(uint8(DecImgCrop2)),  title('Decrypted Image ');

%% Salt and pepper noise attack
fprintf('\nSalt and pepper noise attack ... \n')
figure('Name','Salt and pepper noise attack ...'),
NoiseLevel=0.0005;
EncImgNoise0005=double(imnoise(uint8(EncImg),'salt & pepper',NoiseLevel));
DecImgNoise0005 = Decryption((EncImgNoise0005),Key,H1);
PSNRnosie0005=psnr(uint8(DecImgNoise0005),uint8(PlainImg));
[ssimvalnosie0005] = ssim(uint8(DecImgNoise0005),uint8(PlainImg));
fprintf('Nosiy cipher image Noise Level= %f, PSNR= %f SSIM= %f\n',NoiseLevel,PSNRnosie0005,ssimvalnosie0005);
subplot(2,4,1), imshow(uint8(EncImgNoise0005)),  title('0.0005 Noisy Cipher Image ');
subplot(2,4,5), imshow(uint8(DecImgNoise0005)),  title('0.0005 Noisy Decrypted Image ');
%***********************************************************************************
NoiseLevel=0.005;
EncImgNoise005=double(imnoise(uint8(EncImg),'salt & pepper',NoiseLevel));
DecImgNoise005 = Decryption((EncImgNoise005),Key,H1);
PSNRnosie005=psnr(uint8(DecImgNoise005),uint8(PlainImg));
[ssimvalnosie005] = ssim(uint8(DecImgNoise005),uint8(PlainImg));
fprintf('Nosiy cipher image Noise Level= %f, PSNR= %f SSIM= %f\n',NoiseLevel,PSNRnosie005,ssimvalnosie005);
subplot(2,4,2), imshow(uint8(EncImgNoise005)),  title('0.005 Noisy Cipher Image ');
subplot(2,4,6), imshow(uint8(DecImgNoise005)),  title('0.005 Noisy Decrypted Image ');
%***********************************************************************************
NoiseLevel=0.05;
EncImgNoise05=double(imnoise(uint8(EncImg),'salt & pepper',NoiseLevel));
DecImgNoise05 = Decryption((EncImgNoise05),Key,H1);
PSNRnosie05=psnr(uint8(DecImgNoise05),uint8(PlainImg));
[ssimvalnosie05] = ssim(uint8(DecImgNoise05),uint8(PlainImg));
fprintf('Nosiy cipher image Noise Level= %f, PSNR= %f SSIM= %f\n',NoiseLevel,PSNRnosie05,ssimvalnosie05);
subplot(2,4,3), imshow(uint8(EncImgNoise05)),  title('0.05 Noisy Cipher Image ');
subplot(2,4,7), imshow(uint8(DecImgNoise05)),  title('0.05 Noisy Decrypted Image ');
%***********************************************************************************
NoiseLevel=0.1;
EncImgNoise1=double(imnoise(uint8(EncImg),'salt & pepper',NoiseLevel));
DecImgNoise1 = Decryption((EncImgNoise1),Key,H1);
PSNRnosie1=psnr(uint8(DecImgNoise1),uint8(PlainImg));
[ssimvalnosie1] = ssim(uint8(DecImgNoise1),uint8(PlainImg));
fprintf('Nosiy cipher image Noise Level= %f, PSNR= %f SSIM= %f\n',NoiseLevel,PSNRnosie1,ssimvalnosie1);
subplot(2,4,4), imshow(uint8(EncImgNoise1)),  title('0.1 Noisy Cipher Image ');
subplot(2,4,8), imshow(uint8(DecImgNoise1)),  title('0.1 Noisy Decrypted Image ');
