%------------------------------------------------------------------------------
% Filename: hypgeo_D.m
% 
% To determine confidence interval widths and coverage probabilities of sampling
% from population without replacement - hypergeometric distribution
%
% This code takes population size and sample size to generate sampling
% distributions. Confidence interval widths and coverage probabilities are
% determined for hypergeometric distribution (exact) is
% compared with various approximations: binomial, normal approx and Poisson.
%
% Reference: http://naumaan.com/files/HypGeo-CI-Naumaan_Nayyar.pdf
%
% Author: Naumaan Nayyar
%
% Date: August 8, 2015
%------------------------------------------------------------------------------

close all
clear all
clc

N = 100; %population size
n = 20; %sample size

alpha = 0.05;
alpha1 = 0.025; %left tail probability
alpha2 = alpha - alpha1; %right tail probability

x = 0:n; %number of defectives observed in sample

D = 0:N; %defective size in population
p = zeros(1,N+1); %coverage probability
p2 = zeros(1,N+1); %coverage probability
p3 = zeros(1,N+1); %coverage probability
p4 = zeros(1,N+1); %coverage probability

%% exact confidence interval computation
for i = 1:length(x)
    yu = hygecdf(x(i),N,D,n);
    if(x(i) > 0)
        yl = hygecdf(x(i)-1,N,D,n);
    else
        yl = 0;
    end
    if(yu(end)>alpha1)
        Du(i) = N - (n - x(i));
    else
        Du(i) = find(yu-alpha1 <= 0, 1, 'first')-1; % -1 to match index to zero subscript
        if(Du(i) > 0)
            Du(i) = Du(i) - 1; % -1 to get to needed element
        end
    end
    if(yl(1) < 1 - alpha2)
        Dl(i) = x(i);
    else
        Dl(i) = find(yl-1+alpha2 >= 0, 1, 'last')-1; % -1 to match index to zero subscript
        if(Dl(i) < N)
            Dl(i) = Dl(i) + 1; % +1 to get to needed element
        end
    end
    w(i) = Du(i) - Dl(i);
    p(Dl(i)+1:Du(i)+1) = p(Dl(i)+1:Du(i)+1) + hygepdf(x(i),N,Dl(i):Du(i),n);
end
[x' Dl' Du' w']

%fit((x'-n/2)*2/n,w','poly4') % try to fit a polynomial curve to the width
%as a function of the scaled and centered data

% plot((x-n/2)*2/n,w,'-ob');
% xlabel('Number of defectives observed')
% ylabel('Width of confidence interval of number of defectives in population')
% title('Width of symmetric two-sided confidence interval as a function of number of defectives observed (N=100, n=20)')
% figure,plot(x,Dl,'-og')
% hold on
% plot(x,Du,'-or')
% xlabel('Number of defectives observed')
% ylabel('Lower and upper confidence bound of number of defectives in population')
% title('Two-sided symmetric confidence bounds as a function of number of defectives observed (N=100, n=20)')
% hold off
figure,plot(D(1:end),p(1:end),'-ok');
xlabel('Number of defectives in population');
ylabel('Coverage probability');
title('Coverage probabilities of symmetric two-sided confidence interval as a function of number of defectives observed (N=1000, n=20)');

%% method 2 (naumaan - Normal (hyper) cdf approximation)
contcorr = 0.5; % 0.5 or x/n
c = norminv(alpha1)^2*(N-n)/(N-1);
x = x + contcorr; % continuity correction
Du2 = min(round(N*((2*n*x + c*n) + sqrt(c^2*n^2 + 4*c*n*x.*(n-x)))/(2*(n^2 + c*n))),N);
x = x - contcorr - 0.5; % continuity correction
Dl2 = max(round(N*((2*n*x + c*n) - sqrt(c^2*n^2 + 4*c*n*x.*(n-x)))/(2*(n^2 + c*n))),0);
x = x + 0.5; % restoring x
[x' Dl2' Du2' Du2'-Dl2']
for i = 1:length(x)
    p2(Dl2(i)+1:Du2(i)+1) = p2(Dl2(i)+1:Du2(i)+1) + hygepdf(x(i),N,Dl2(i):Du2(i),n);
end
figure,plot(D(2:end-1),p2(2:end-1),'-ok');
xlabel('Number of defectives in population');
ylabel('Coverage probability');
title('Coverage probabilities of symmetric two-sided confidence interval [Normal (hyper)] (N=100, n=20)');

%% method 3 (normal error approximation Wald-type)
fpc = sqrt((N-n)/(N-1));
x = x/n; % proportion conversion
Dl3 = max(round(N*(x + fpc*norminv(alpha1)*sqrt(x.*(1-x)/n))),0);
Du3 = min(round(N*(x - fpc*norminv(alpha1)*sqrt(x.*(1-x)/n))),N);
x = x*n; % restoring x
[x' Dl3' Du3' Du3'-Dl3']
for i = 1:length(x)
    p3(Dl3(i)+1:Du3(i)+1) = p3(Dl3(i)+1:Du3(i)+1) + hygepdf(x(i),N,Dl3(i):Du3(i),n);
end
figure,plot(D(1:end),p3(1:end),'-ok');
xlabel('Number of defectives in population');
ylabel('Coverage probability');
title('Coverage probabilities of symmetric two-sided confidence interval [Normal (Wald-type)] (N=100, n=20)');


%% method 4 - clopper-pearson binomial interval (exact computation using inverse of beta cdf)
Du4 = floor(min(N*betainv(1-alpha1,x+1,n-x),N));
Dl4 = ceil(max(0,N*betainv(alpha2,x,n-x+1)));
[x' Dl4' Du4' Du4'-Dl4']
for i = 1:length(x)
    p4(Dl4(i)+1:Du4(i)+1) = p4(Dl4(i)+1:Du4(i)+1) + hygepdf(x(i),N,Dl4(i):Du4(i),n);
end
figure,plot(D(1:end),p4(1:end),'-ok');
xlabel('Number of defectives in population');
ylabel('Coverage probability');
title('Coverage probabilities of symmetric two-sided confidence interval [Binomial] (N=100, n=20)');

figure,plot(x,w,'-ob',x,Du2-Dl2,'-xr',x,Du3-Dl3,'-*k',x,Du4-Dl4,'-+g');
xlabel('Number of defectives observed')
ylabel('Width of confidence interval of number of defectives in population')
title('Width of symmetric two-sided confidence interval as a function of number of defectives observed (N=100, n=20)')
legend('exact','Normal (hyper)','Normal (Wald-type)','binomial')