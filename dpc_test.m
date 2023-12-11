%% khanna dpc test
clear
close all
load('N:\WissPers\Users\Chan\PhD\Signale\NFT DP\Keysight_dp_0+0-3i0+0-6i_num_0_2ev_4psk_0_45_phaseshift.mat')
error_max=3;
Length_Signal=length(ref_signal_x);

%% amplitude fluctuations/noise
amp_rand=randi([9000 10000], Length_Signal,1)/10000;
% amp_rand=1;
signal_rx_x=ref_signal_x.*amp_rand.';
amp_rand=randi([9000 10000], Length_Signal,1)/10000;
% amp_rand=1;
signal_rx_y=ref_signal_y.*amp_rand.';

%% IQ Imbalance % massive destruction of EVs and b-coeffs
I_imbalance_x=0.98;
Q_imbalance_x=1.2;
I_imbalance_y=1.00;
Q_imbalance_y=1.1;

signal_rx_x=real(signal_rx_x)*I_imbalance_x+1i*imag(signal_rx_x)*Q_imbalance_x;
signal_rx_y=real(signal_rx_y)*I_imbalance_y+1i*imag(signal_rx_y)*Q_imbalance_y;

% phase rotation % mainly influences the constellation phase. Can also seperate constellation points into clouds!
phase_rot = deg2rad(25);
jitter_offset=5;
phase_jitter=deg2rad(-jitter_offset+(2*jitter_offset*rand(1,length(signal_rx_x))));
signal_rx_x=signal_rx_x.*exp(1i*phase_rot+phase_jitter);

phase_rot = deg2rad(64);
jitter_offset=7;
phase_jitter=deg2rad(-jitter_offset+(2*jitter_offset*rand(1,length(signal_rx_y))));
signal_rx_y=signal_rx_y.*exp(1i*phase_rot+phase_jitter);

% nonlinearity with sine
signal_rx_x=sin(signal_rx_x);
signal_rx_y=sin(signal_rx_y);


%% %%%%%%%%% NFT
t=param.t;
siminftnft=1;
 UpsamplingFactor=5;
 Num_Symbols=param.num_sym-param.n_noinfo;
 param.nft_method=1;
 EVs=param.EVs;
 t_norm=param.t_norm;

 signal_rx_x=signal_rx_x(1:Num_Symbols*param.num_samples);
 signal_rx_y=signal_rx_y(1:Num_Symbols*param.num_samples);

    signal_nft_x=resample(signal_rx_x,UpsamplingFactor,1);
    signal_nft_y=resample(signal_rx_y,UpsamplingFactor,1);

    signal_nft_x=reshape(signal_nft_x,1,UpsamplingFactor*N_Samples,Num_Symbols);
    signal_nft_y=reshape(signal_nft_y,1,UpsamplingFactor*N_Samples,Num_Symbols);
    
    param.t=linspace(t(1),t(end),param.num_samples*UpsamplingFactor);
    param.num_samples=param.num_samples*UpsamplingFactor;
    
    [b1temp,b2temp,~,~,~,~,~,~,eigenvaluestemp, param, out.nft_time]=NFT_ken_dp(signal_nft_x, signal_nft_y,param,siminftnft);

        b1=b1temp;
        b2=b2temp;
        eigenvalues=eigenvaluestemp;
        
    % end
        figure;
        set(gcf, 'Units', 'pixel', 'Position', [1920, 50, 1800 900])
        tiledlayout(2,3)

        nexttile
        hold on
        for j=1:length(EVs)
            plot(eigenvalues(:,j),'.')
        end
        axis([-.5 .5 0 1])
        pbaspect([1 1 1])
        title('Eigenwerte')
        grid minor
        hold off

        nexttile
        hold on
        for j=1:length(EVs)
            plot(b1(:,j),'.')
        end
        pbaspect([1 1 1])
        axis([-1 1 -1 1]*max(param.abs_q_disc))
        title('b1')
        grid minor
        hold off
        
        nexttile
        hold on
        for j=1:length(EVs)
            plot(b2(:,j),'.')
        end
        pbaspect([1 1 1])
        axis([-1 1 -1 1]*max(param.abs_q_disc))
        title('b2')
        grid minor
        hold off
        
        nexttile % empty
        pbaspect([1 1 1])

        ylim_border=max(max([abs(signal_nft_x) abs(signal_nft_y)]));
        t_resampled=linspace(-t_norm,t_norm,N_Samples*UpsamplingFactor);
        nexttile
        hold on
        for j=1:floor(Num_Symbols/2)
        plot(t_resampled,abs(signal_nft_x(1,:,j)))
        end
        pbaspect([1 1 1])
        title('Abs. val. of solitons x-pol')
        xlim([-t_norm t_norm])
        ylim([0 ylim_border])
        grid minor

        nexttile
        hold on
        for j=1:floor(Num_Symbols/2)
        plot(t_resampled,abs(signal_nft_y(1,:,j)))
        end
        pbaspect([1 1 1])
        title('Abs. val. of solitons y-pol')
        xlim([-t_norm t_norm])
        ylim([0 ylim_border])
        grid minor


X=[real(ref_signal_x(1:Num_Symbols*param.num_samples/UpsamplingFactor)); imag(ref_signal_x(1:Num_Symbols*param.num_samples/UpsamplingFactor)); real(ref_signal_y(1:Num_Symbols*param.num_samples/UpsamplingFactor)); imag(ref_signal_y(1:Num_Symbols*param.num_samples/UpsamplingFactor))];
%% combine rx signal
Y=[real(signal_rx_x); imag(signal_rx_x); real(signal_rx_y); imag(signal_rx_y)];

%% parameters
P=3; % nonlinearities
mem=5; % memory (should be odd number)
N_use=120*10;
N=N_use+mem-1;

clear h H h_use h_use_str X_temp
% h=zeros(mem*P,4,dpc_trial_max+1); % max_iter+1 since we have the first entry as initialization

h=zeros(mem*P,4,error_max+1); % max_iter+1 since we have the first entry as initialization
h(1,1:4,:)=1;
err_real=zeros(error_max+1,1);

for zz=1:4
    clear x y y_new z z_plot err h_use h_use_str_temp min_err err_real err_real_str z_full y_full
    y=Y(zz,1:N);
    x=X(zz,1:N);

 figure
     hold on
     xlim([1 1000])

    for errortrial=1:error_max
        %% misc

        %% DPC p.53 fig. 4.3
        h_use=reshape(squeeze(h(:,zz,errortrial)),mem,P);

        z=mpm_nonlin(y, h_use, N, P, mem);

        z=[zeros(1, floor(mem/2)) z zeros(1, floor(mem/2))];

plot(z)

        % output AFTER channel/ optical transmission system (the output must be feeded through the channel again after compensation)
        [h(:,zz,errortrial+1), err_real(errortrial+1,:)]=ila_nonlin(x, y, h_use, N_use, P ,mem);

        %% update y
        y=z;
        disp(sprintf('EVM between z and z_hat for channel %d run %d', zz, errortrial))
        disp(err_real(errortrial+1,:));

        err_real_str=cell(zz+2,1);
        err_real_str{1}='rx';
        err_real_str{2}='tx';
        for i=3:size(err_real,1)+2
            err_real_str{i}=num2str(round(err_real(i-2,:),3));
        end
    end

    %% get h-coefficients from lowest error run
    err_real=err_real(2:end);
    [~,min_err]=min(err_real);
    %         h=h(:,2:end);
    h_use=h(:,zz,min_err+1);

    h_use_str_temp=[];
    for i=1:size(h_use,1)
        h_use_str_temp=[h_use_str_temp ' ' num2str(h_use(i))];
    end
    %
    h_use_str{zz,:}=h_use_str_temp;

    %% full signal
    % only use linear tap of 3rd and 5th order as in (5.4) p.89
    %     h_use_new=reshape(squeeze(h(:,zz,error_max+1)),mem,P);

    %     h_use_new(1,1)=1;
    h_use_new=reshape(h_use,mem,P);

    % since we use memorial polynomial model, we need to zeropad the
    % original signal with mem/2 after the signal
    x_temp=[X(zz,:) zeros(1,2*floor(mem/2))];

    signal_DPCed_temp=mpm_nonlin(x_temp, h_use_new, Length_Signal+(mem-1), P, mem);
%     signal_DPCed(zz,:)=rescale(signal_DPCed_temp,-1,1);
    signal_DPCed(zz,:)=signal_DPCed_temp;

    H(zz,:)=h_use;
end

%% plot inputs/outputs
a=figure;
a.Position = [1930 20 1900 950];
tiledlayout(2,2)
for zz=1:4
    nexttile
    hold on
    plot(X(zz,:),'k','LineWidth',1.25)
    plot(Y(zz,:),'r','LineWidth',1.25)
    plot(signal_DPCed(zz,:),'-.','LineWidth',1.25)
    title(sprintf('TX, RX, and optimal output y with h-coefficients %s',h_use_str{zz,:}))
    legend('TX','RX', 'precompensated')
    ylim([-1 1]*1.5)
    xlabel('Normed amp.')
    ylabel('Samples')
    xlim([1 3000]);
    grid on
end

%% save new input data
% ref_signal_x_dpc=X_DPCed(1,:)+1i*circshift(X_DPCed(2,:),-1);
% ref_signal_y_dpc=circshift(X_DPCed(3,:),2)+1i*circshift(X_DPCed(4,:),1);

signal_x=signal_DPCed(1,:)+1i*signal_DPCed(2,:);
signal_y=signal_DPCed(3,:)+1i*signal_DPCed(4,:);

%% %%%%%%%%% NFT
t=param.t;
siminftnft=1;
 UpsamplingFactor=5;
 Num_Symbols=param.num_sym-param.n_noinfo;
 param.nft_method=1;
 EVs=param.EVs;

    signal_nft_x=resample(signal_x(1:Num_Symbols*param.num_samples),UpsamplingFactor,1);
    signal_nft_y=resample(signal_y(1:Num_Symbols*param.num_samples),UpsamplingFactor,1);

    signal_nft_x=reshape(signal_nft_x,1,UpsamplingFactor*N_Samples,Num_Symbols);
    signal_nft_y=reshape(signal_nft_y,1,UpsamplingFactor*N_Samples,Num_Symbols);
    
    param.t=linspace(t(1),t(end),param.num_samples*UpsamplingFactor);
    param.num_samples=param.num_samples*UpsamplingFactor;
    
    [b1temp,b2temp,abs_q1,abs_q2,q_cont1,q_cont2,q_cont1_mean_energy,q_cont2_mean_energy,eigenvaluestemp, param, out.nft_time]=NFT_ken_dp(signal_nft_x, signal_nft_y,param,siminftnft);

        b1=b1temp;
        b2=b2temp;
        eigenvalues=eigenvaluestemp;
        
    % end
        figure;
        set(gcf, 'Units', 'pixel', 'Position', [1920, 50, 1800 900])
        tiledlayout(2,3)

        nexttile
        hold on
        for j=1:length(EVs)
            plot(eigenvalues(:,j),'.')
        end
        axis([-.5 .5 0 1])
        pbaspect([1 1 1])
        title('Eigenwerte')
        grid minor
        hold off

        nexttile
        hold on
        for j=1:length(EVs)
            plot(b1(:,j),'.')
        end
        pbaspect([1 1 1])
        axis([-1 1 -1 1]*max(param.abs_q_disc))
        title('b1')
        grid minor
        hold off
        
        nexttile
        hold on
        for j=1:length(EVs)
            plot(b2(:,j),'.')
        end
        pbaspect([1 1 1])
        axis([-1 1 -1 1]*max(param.abs_q_disc))
        title('b2')
        grid minor
        hold off
        
        nexttile % empty
        pbaspect([1 1 1])

        ylim_border=max(max([abs(signal_nft_x) abs(signal_nft_y)]));
        t_resampled=linspace(-t_norm,t_norm,N_Samples*UpsamplingFactor);
        nexttile
        hold on
        for j=1:floor(Num_Symbols/2)
        plot(t_resampled,abs(signal_nft_x(1,:,j)))
        end
        pbaspect([1 1 1])
        title('Abs. val. of solitons x-pol')
        xlim([-t_norm t_norm])
        ylim([0 ylim_border])
        grid minor

        nexttile
        hold on
        for j=1:floor(Num_Symbols/2)
        plot(t_resampled,abs(signal_nft_y(1,:,j)))
        end
        pbaspect([1 1 1])
        title('Abs. val. of solitons y-pol')
        xlim([-t_norm t_norm])
        ylim([0 ylim_border])
        grid minor
