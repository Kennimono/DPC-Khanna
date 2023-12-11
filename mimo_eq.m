function [Y,C]=mimo_eq(X,Y_hat)
%% MIMO Equalizer from Khanna diss p.69 chapter 4.5.1
% based on model p.66 fig.4.11 and (4.40) 
% after upsampling

% x      original tx samples of each channel                   size N x 1
% X      original tx 4-D signal                                size N x 4
% y      estimated output samples                              size N x 1
% Y      estimated de-multiplexed 4-D signal, input for ILA    size N x 4
% Y_hat  mixed channels 4-D signal                             size N x 4
% c      channel coefficients of DP-MZM                        size 1x1
% C      channel matrix                                        size 4x4
% C_hat  estimated channel matrix                              size 4x4
% c      coefficients with minimalized errors                  size 1x1
% r      cross correlation vector between y_hat(i) and x(k)    size 4x4
% tau    delays between k-th txed and i-th rx-ed channel       size 1x1
% Tau    delay matrix with tau                                 size 4x4
% T0     Sampleperiod of ADC/Scope                             size 1x1
%
% n      number of samples
% N      max number of samples
% k      k-th transmitted channel
% i      i-th received channel

%% initialization
N=size(Y_hat,1);

r=zeros(4,4,N);
Tau=zeros(4,4,N);
C=zeros(4,4,N);
Y=zeros(4,4,N);

%% after rx, Y_hat consists of 4 channel signals in the form p.69 (4.47)
% y_hat(i)=sum k=1:4 c(i,k)*y_k(n*T0-tau(i,k))+w_i(n*T0)

%% coarse delay estimation and compensation based on (4.34-4.36)
% for i=1:4
%     for k=1:4
% %         r_temp=fft((fliplr(Y_hat(k,:)))).*fft(conj(X(k,:))); % original deifniton ohne ifft
%         r_temp=ifft(fft((fliplr(Y_hat(k,:)))).*fft(conj(X(k,:))));       
% %         [~,max_ind(i,k)]=max(r(i,k,:));
% [~,Tau(i,k)]=max(r_temp);
%         
% %         Tau(i,k)=r_temp(max_ind(i,k));
%     end
% end

%% finding original Y with p.70 (4.51-4.53)
C_hat=X*X.'\Y_hat*X.';

% C_hat=C_hat./max(abs(C_hat));

Y=C_hat\Y_hat; % using pseudo inverse

% for i=1:4
%     Y(i,:)=Y(i,:)/max(abs(Y(i,:)));
% end

% figure
% tiledlayout(2,4)
% nexttile
% plot(Y_hat(1,:))
% ylim([-1 1]*1.1)
% nexttile
% plot(Y_hat(2,:))
% ylim([-1 1]*1.1)
% nexttile
% plot(Y_hat(3,:))
% ylim([-1 1]*1.1)
% nexttile
% plot(Y_hat(4,:))
% ylim([-1 1]*1.1)
% nexttile
% plot(Y(1,:))
% ylim([-1 1]*1.1)
% nexttile
% plot(Y(2,:))
% ylim([-1 1]*1.1)
% nexttile
% plot(Y(3,:))
% ylim([-1 1]*1.1)
% nexttile
% plot(Y(4,:))
% ylim([-1 1]*1.1)