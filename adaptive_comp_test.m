%% khanna test with 1 channel qam
clear
close all
warning('off','all')

M = 64; % M-QAM
phaseshift_vec = deg2rad([0 0]) ; % phaseshift of qam
% x = (0:M-1)';
Length_Signal=96000;
span = 1; % Filter span in symbols (signal delayed by n symbols)
SamplerateTX = 120e9; % How many Samples per second
T_Sample_tx=1/SamplerateTX; % Sample Duration
Baudrate = 1e9; % Samplerate/n, n = actual Samplenumber per symbol
T_Sym=1/Baudrate;
n_samples_tx = SamplerateTX/Baudrate; % Number of samples
symbols = floor(Length_Signal/n_samples_tx);  % number of synbols (later last n symbols will be set to 0 for correlation)

racos_alpha=1; % roll off factor for rasied cosine
FreqVec_tx=(-Length_Signal/2:Length_Signal/2-1)/(Length_Signal*T_Sample_tx);

%% generate data
    rng(1337)
    data_tx_x = randi([0 M-1],symbols,1);

sig_temp_x = qammod(data_tx_x,M, 'bin', 'UnitAveragePower',true).*exp(1i*phaseshift_vec(1));

%% Pulse shape signal and generate signal
txfilter = comm.RaisedCosineTransmitFilter('RolloffFactor',racos_alpha, ...
    'FilterSpanInSymbols',span,'OutputSamplesPerSymbol',n_samples_tx);

% generate time domain signal and add zeros for sequence correlation at rx
signal_x = txfilter(sig_temp_x).';
signal_x=rescale(real(signal_x),-1,1)+1i*rescale(imag(signal_x),-1,1);

%% demod to check
rxfilter = comm.RaisedCosineReceiveFilter('Shape', 'Normal', 'RolloffFactor',racos_alpha, ...
    'FilterSpanInSymbols',span,'InputSamplesPerSymbol',n_samples_tx, ...
    'DecimationFactor',n_samples_tx);
sig_rx_x=rxfilter((signal_x./sqrt(M)).');

%% parameters
P=3; % nonlinearity order
% L=[1 0 3 0 5];
L=[1 0 3];
mem=5; % memory (odd number)
N=48000;
% N=N_use+2*floor(mem/2);

%% misc
max_iter=20; % maximal iterations

% N_dpc=N-2*M;
%% inputs follow figure 4.3 p. 53
% input to system

x=real(signal_x(1:N+mem-1));
h=zeros(mem*P,1); % max_iter+1 since we have the first entry as initialization
h(1)=1;
    
%% DPC p.53 fig. 4.3 (indirect learning architecture)
% while err > 1 && max_iter~=0
for q=1:max_iter
    clear y y_new z z_plot err h_use_str min_err err_real_str z_full y_full % clear after every iteration
    if q==1
h_use=h;
    end

    % we dont needd mpm at all, just filter with h coefficient!
    z=filter(h_use,1,x);
        % if q==1
        % z=mpm_nonlin(x, h_use, N, mem, P, L);
        % else
        % z=mpm_nonlin(y_comp, h_use, N, mem, P, L);
        % end
%         z=[zeros(1, floor(mem/2)) z zeros(1, floor(mem/2))];

% output AFTER channel/ optical transmission system (the output must be feeded through the channel again after compensation)
% rng(1337)
%% awgn channel
SNR=10;
y=awgn(z,SNR);
y=asin(y);
y=real(y);

% normalize signals
z=z/max(abs(z));
y=y*sqrt(sum(z.^2)/sum(y.^2));

[h(:,q+1), err_real(q,:)]=ila_nonlin(x, y, N, P ,mem);
% [h(:,q+1), err_real(q,:)]=indirect_learning_architecture(x, z, y, h(:,q), N_use, P ,M);

% max_iter=max_iter-1;
y_comp=mpm_nonlin(y, reshape(h(:,q+1),mem,P), N, mem, P, L);

%% we use the error between x and y_new, since this is the real metric for the coefficient.
% err_real(q,:)=norm(x-y_new(:,q).');
% q=q+1;
y_new(q,:)=y_comp;
h_use=h(:,q+1);
end
disp('EVM between x and y_new')
disp(err_real);

err_real_str=cell(q+2,1);
err_real_str{1}='rx';
err_real_str{2}='tx';
for i=3:size(err_real,1)+2
err_real_str{i}=num2str(round(err_real(i-2,:),3));
end

%% get h-coefficients from lowest error run
[~,min_err]=min(err_real);
h_use=h(:,min_err+1);

h_use_str=[];
for i=1:size(h_use,1)
h_use_str=[h_use_str ' ' num2str(round(h_use(i,:),3))];
end

%% full compensated signal
z_full=conv(real(signal_x),h_use,'same');

%% plot inputs/outputs
a=figure;
tiledlayout(3,1)
a.Position = [1930 20 800 900];
% nexttile
% hold on
% for i=1:max_iter
% plot(z_plot(:,i))
% end
% title('Compensated input data z with corresponding EVMs')
% legend(err_real_str{3:end})
% ylim([-1 1]*1.5)
% xlabel('Normed amp.')
% ylabel('Samples')

nexttile
plot(y,'LineWidth',2)
hold on
plot(x,'.-k','LineWidth',2)
for i=1:max_iter
plot(y_new(i,:))
end
title('Original TX and all compensated output data y with corresponding EVMs')
legend (err_real_str)
% ylim([-1 1]*1.5)
xlabel('Normed amp.')
ylabel('Samples')
xlim([1 1000])

nexttile
hold on
plot(x,'.-k','LineWidth',2)
plot(y_new(min_err,:),'LineWidth',2)
title(sprintf('Original x and optimal output y with h-coefficients %s',h_use_str))
legend('orig', 'compensated')
% ylim([-1 1]*1.5)
xlabel('Normed amp.')
ylabel('Samples')
xlim([1 1000])

nexttile
hold on
plot(real(signal_x),'.-k','LineWidth',2)
plot(z_full,'LineWidth',2)
title(sprintf('Full-length x and optimal y with h-coefficients %s',h_use_str))
legend('orig', 'compensated')
% ylim([-1 1]*1.5)
xlabel('Normed amp.')
ylabel('Samples')
xlim([1 1000])
