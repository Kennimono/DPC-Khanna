function y=mpm_nonlin_causal(x, h, N, P, M)
%% memory polynomial model with 3rd or 5th order
% from khanan diss p. 89 (5.4)
% can be modelled as time invariant non-lienar system with memory (traditionally modelled as volterra series)

% n current number of sample
% N max number of samples
% x input sample                     size 1 x N
% y output sample                    size 1 x N
% h DPC coefficients                 size 1 x (M+1)*P
% m current memory tap
% M max memory size
% p p-th order impulse response (IR)
% P max IR order

%% initialize stuff
y=zeros(1,N-M+1);

%% actual model
% in matlab you cant have negative index. so we just shift the
% representation of (4.6) M/2 to the right to get the non causal structure
l=1;
for n=1+floor(M/2):N-floor(M/2)
    y_temp=0;
    for m=-floor(M/2):1:floor(M/2)
        y_temp=y_temp+h(floor(M/2)+m+1,1)*x(n+m);
    end
    if P==3
        y(l)=y_temp+x(n).^3.*h(ceil(M/2),3);
    elseif P ==5
        y(l)=y_temp+x(n).^3.*h(ceil(M/2),3)+x(n).^5.*h(ceil(M/2),5);
    else % linear
        y(l)=y_temp;
    end
    l=l+1;
end

% y=y/max(abs(y));
% figure
% hold on
% plot(x)
% plot(y)