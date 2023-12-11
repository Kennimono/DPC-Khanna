function [h, err, z_hat]=ila_nonlin_3_5_linear(x, z, y, h, N_use, P, M)
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
% y=y*sqrt(sum(z.^2)/sum(y.^2));
y=y*sqrt(sum(x.^2)/sum(y.^2));
% z=x;
B_m=ones(N_use,M);
% z_hat=zeros(N_use,1);

%% step 2: memorial polynomial model with M memory taps and order P
B=[];
for p=1:P
    l=1;
    for m=floor(M/2):-1:-floor(M/2)
        B_m(:,l)=[y(floor(M/2)-m+1:N_use+floor(M/2)-m)].^p;
        l=l+1;
    end
    B=[B B_m];
end

%% step 3: find coefficients h, that the error err is minimized

% normal equation; eliminate alldependent columns in B beforehand
if isinf(inv(B.'*B))
    error('Linear dependant rows/columns in inv(conj(B)*B) -.-')
end

h_hat=inv(B.'*B)*B.'*x(1+floor(M/2):end-floor(M/2));   % be careful to do matrix multiplication with B, not elementwise!
% h_inverse=inv(B.'*B)*B.'*x(1+floor(M/2):end-floor(M/2));

% calculate error
err=norm(x(1+floor(M/2):end-floor(M/2))-B*h_hat);
% err=min(norm(z-z_hat));

%% step 4 obtain h with following conditions. h cxan be then convoluted with x to obtain predistorted signal z
% x->y, z->z_hat, e->0
h=h_hat;