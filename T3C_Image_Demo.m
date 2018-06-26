% addpath
currentFolder = pwd;
addpath(genpath(currentFolder));
%% Data preparation
image=imread('lena.bmp');
image=imresize(image,[256 256]);
T=double(image)/255; % regularization
S=size(T); N=numel(S);
mr=0.5; % missing rate
W=gen_W(S,mr);
X=T.*W;
r=[150,3];

%% Run two algorithms

% TTWOPT- quick convergence
[X_hat_wopt,G_wopt]=T3C(X,W,'Alg','TTWOPT','Rank',r,'MaxIter', 1e3,'Tol', 1e-4);

% TTSGD - able to deal large-scale data and low complexity
[X_hat_sgd,G_sgd]=T3C(X,W,'Alg','TTSGD','Rank',r,'MaxIter', 5e5,'Tol', 1e-4);

%% Performance evaluation 

RSE_wopt=RSE_fun(T,X_hat_wopt,W);
RSE_sgd=RSE_fun(T,X_hat_sgd,W);
PSNR_wopt = PSNR_RGB(double(255*X_hat_wopt),double(image));
SSIM_wopt = ssim_index(rgb2gray(uint8(255*X_hat_wopt)),rgb2gray(uint8(image)));
PSNR_sgd = PSNR_RGB(double(255*X_hat_sgd),double(image));
SSIM_sgd = ssim_index(rgb2gray(uint8(255*X_hat_sgd)),rgb2gray(uint8(image)));

%% Report results
figure('position', [250, 400, 900, 300]);
subplot(1,4,1);  imshow(image);title('Original');set(gca,'Fontsize',20);
subplot(1,4,2);  imshow(uint8(255*X));title('Observed');set(gca,'Fontsize',20);
subplot(1,4,3);  imshow(uint8(255*X_hat_wopt));title('TTWOPT');set(gca,'Fontsize',20);
subplot(1,4,4);  imshow(uint8(255*X_hat_sgd));title('TTSGD');set(gca,'Fontsize',20);
fprintf('-----Evaluation reports-----\n\n')
fprintf('Tensor size');disp(S);
fprintf('\n Missing rate=%g\n\n',mr);
fprintf('       \t   RSE \t PSNR \t SSIM\n TTWOPT   %.3f\t %.2f \t %.3f \t \n TTSGD    %.3f  %.2f \t %.3f\t\n ',RSE_wopt(1),PSNR_wopt,SSIM_wopt,RSE_sgd(1),PSNR_sgd,SSIM_sgd);
fprintf('\n----------------------------\n')