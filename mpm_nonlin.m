function y=mpm_nonlin(x, h, N, M, P, L)
%% memory polynomial model with 3rd or 5th order
%  23-11-21 ver.
% from khanan diss p. 89 (5.4)
% can be modelled as non-causal time invariant non-linear system with memory (traditionally modelled as volterra series)

% n current number of sample
% N max number of samples
% x input sample                     size 1 x N
% y output sample                    size 1 x N+M
% h DPC coefficients                 size 1 x (M+1)*P
% m current memory tap
% M max memory size
% p p-th order impulse response (IR)
% P max IR order
% L weight of order P

%% initialize stuff
y=zeros(1,N+M-1);

%% actual model
% in matlab you cant have negative index. so we just shift the
% representation of (4.6) M/2 to the right to get the non causal structure;
for p=1:P
y=y+conv(h(:,p),x(1:N).^L(p));
end

% y=y/max(abs(y));
% figure
% hold on
% plot(x)
% plot(y)