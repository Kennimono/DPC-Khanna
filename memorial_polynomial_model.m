function y=memorial_polynomial_model(x, N, M, P)
%% memory polynomial model
% from khanna diss p. 51 ff chapter 4.1
% linear model for 3dB-bandwidth effect of transmitter and I/Q-skew (DAC, amp, DP-MZM)
% can bemodelled as time invariant non-lienar system with memory (traditionally modelled as volterra series)

% n number of samlpes
% N max number of samples
% x input sample                     size 1 x N
% y output sample                    size 1 x N
% h volterra kernels, (coefficients) size 1 x (M+1)*P
% m memory tap
% M max memory size
% p p-th order impulse response (IR)
% P max IR order

%% initialize stuff
% TODO

h=ones(M,P);
y=x;

%% actual model
% in matlab you cant have negative index. so we just shift the
% representation of (4.6) M/2 to the right to get the non causal structure
for p=1:P
    for n=round(M/2):N-round(M/2)
        y_temp=y(n);
        for m=1:M
            y_temp=y_temp+h(m,p)*x(n-(M/2)+m).^p;
        end
        y(n)=y_temp;
    end  
end