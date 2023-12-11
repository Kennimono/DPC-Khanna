function [h, err, z_hat]=ila_nonlin_old(x, z, y, h, N_use, P, M)
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

if size(z,1) <= size(z,2)
z=z.';
end

if size(y,1) <= size(y,2)
y=y.';
end

if size(h,1) <= size(h,2)
h=h.';
end

H=reshape(h,M,P).';

%% step 1: collect N samples of x and y and norm
x=x/max(abs(x));
z=z/max(abs(z));
y=y*sqrt(sum(z.^2)/sum(y.^2));


% z=x;
B_m=ones(N_use,M);
z_hat=zeros(N_use,1);

% check h dimension
% if size(h,1)~=M && size(h,2)~=P
%     h=reshape(h,M,P);
% end

%% step 2: memorial polynomial model with M memory taps and order P
% eq 4.11 (dont need to use)
% z_hat is output of training block
% for n=1:N_use
%     l=1;
%     for m=floor(M/2):-1:-floor(M/2)
%         z_hat(n)=z_hat(n)+h(l)*y(n+floor(M/2)-m);
%         l=l+1;
%     end
%     
%     if P == 3
%         z_hat(n)=z_hat(n)+H(3,floor(M/2))*z(n+floor(M/2)).^3;
%     elseif P == 5
%         z_hat(n)=z_hat(n)+H(3,floor(M/2))*z(n+floor(M/2)).^3+H(5,floor(M/2))*z(n+floor(M/2)).^5;
%     elseif mod(P,2)==0 && P > 1
%         for zz=1:P
%         z_hat(n)=z_hat(n)+H(zz,floor(M/2))*z(n+floor(M/2)).^zz;
%         end
%     end
% end

B=[];
for p=1:P
    for m=1:M
        %         B_m(:,m)=circshift(y.'.^p,m-1);
        B_m(:,m)=[y(m:N_use-1+m)].^p;
    end
    B=[B B_m];
end

z_hat=B.*reshape(h,1,M*P);
z_hat=z_hat(:,floor(M/2)+1); % no delay

%% step 3: find coefficients h, that the error err is minimized

% normal equation; eliminate alldependent columns in B beforehand
% B_eliminated=B(:,1);
% h_hat=(inv((B_eliminated.'.*B_eliminated)).*B_eliminated.'*x.').';
if isinf(inv(B.'*B))
    error('Linear dependant rows/columns in inv(conj(B)*B) -.-')
end

 h_hat=inv(B.'*B)*B.'*z(1+floor(M/2):end-floor(M/2));   % be careful to do matrix multiplication with B, not elementwise!
% end

% calculate error
err=norm(z(1+floor(M/2):end-floor(M/2))-B*h_hat);
% err=min(norm(z-z_hat));

%% step 4 obtain h with following conditions. h cxan be then convoluted with x to obtain predistorted signal z
% x->y, z->z_hat, e->0
h=h_hat;