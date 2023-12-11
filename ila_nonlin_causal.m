function [h, err]=ila_nonlin_causal(x, z, y, h, N_use, P, M)
% warning('off','all')

%% indirect learning architecture
% from khanna diss p 52 ff. chapter 4.2
% I dont use (M+1)*P dimension because indexing starts at 1 in matlab

% x      original tx samples                                   size N x 1
% y      output samples after channel/system                   size N x 1
% z      predistorted samples into channel/system              size N x 1
% z_hat  predistorted samples from training; z_hat=B*h         size N x 1
% err    error vector                                          size N x 1
% h      volterra kernels, (coefficients)                      size M*P x 1
% h_hat  coefficients with minimalized errors                  size M*P x 1
% H      transfer function model for coefficients h
% B      H=B*h_hat, B with time shifted copies of y    size N x M*P
%
% n      number of samples
% N      max number of samples
% N_use  used number of samples
% m      memory tap
% M      max memory size
% p      p-th order impulse response (IR)
% N_dpc  extracted sequence with +-M distance to beginnign and end for memory
% P      max IR order

%% initial check: dimensions correction
if size(x,1) <= size(x,2)
x=x.';
end

% if size(z,1) <= size(z,2)
% z=z.';
% end

if size(y,1) <= size(y,2)
y=y.';
end

if size(h,1) <= size(h,2)
h=h.';
end

H=reshape(h,M,P).';

%% step 1: collect N samples of x and y and norm
x=x/max(abs(x));
% z=z/max(abs(z));
y=y*sqrt(sum(x.^2)/sum(y.^2));
% x_use=x(1+M:end-M);
% y_use=y(1+M:end-M);
% figure
% hold on
% plot(x)
% plot(y)

% z=x;
B_m=ones(N_use,M);
% z_hat=zeros(N_use,1);

%% step 2: memorial polynomial model with M memory taps and order P
B=[];
for p=1:P
    for m=1:M
        %         B_m(:,m)=circshift(y.'.^p,m-1);
        B_m(:,m)=y(m:N_use+m-1).^p;
    end
    B=[B B_m];
end

% z_hat=B.*reshape(h,1,M*P);
% z_hat=z_hat(:,floor(M/2)+1); % no delay

%% step 3: find coefficients h, that the error err is minimized

% normal equation; eliminate alldependent columns in B beforehand
% B_eliminated=B(:,1);<
% h_hat=(inv((B_eliminated.'.*B_eliminated)).*B_eliminated.'*x.').';
if isinf(inv(B.'*B))
    error('Linear dependant rows/columns in inv(conj(B)*B) -.-')
end

h_hat=inv(B.'*B)*B.'*x(1:N_use);   % be careful to do matrix multiplication with B, not elementwise!

% calculate error
err=norm(x(1:N_use)-B*h_hat);
% err=min(norm(z-z_hat));

%% step 4 obtain h with following conditions. h cxan be then convoluted with x to obtain predistorted signal z
% x->y, z->z_hat, e->0
h=h_hat;