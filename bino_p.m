%------------------------------------------------------------------------------
% Filename: bino_p.m
% 
% To determine confidence interval widths and coverage probabilities of sampling
% from population with replacement - binomial distribution
%
% This code takes population size and sample size to generate sampling
% distributions. Confidence interval widths and coverage probabilities are
% determined for binomial distribution.
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

N = 1000; %population size
n = 20; %sample size

alpha = 0.05;
alpha1 = 0.025; %left tail probability
alpha2 = alpha - alpha1; %right tail probability

x = 0:n; %number of defectives observed in sample

D = 0:N; %defective size in population
p = zeros(1,N+1); %coverage probability

%% method 1 - exact computation
for i = 1:length(x)
    yu = binocdf(x(i),n,D/N);
    if(x(i) > 0)
        yl = binocdf(x(i)-1,n,D/N);
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
    p(Dl(i)+1:Du(i)+1) = p(Dl(i)+1:Du(i)+1) + binopdf(x(i),n,(Dl(i):Du(i))/N);
end
[x' Dl' Du' w']

%fit((x'-n/2)*2/n,w','poly4') % try to fit a polynomial curve to the width
%as a function of the scaled and centered data

% plot(x,w,'-ob');
% xlabel('Number of defectives observed')
% ylabel('Width of confidence interval of number of defectives in population')
% title('Width of symmetric two-sided confidence interval as a function of number of defectives observed (N=1000, n=50)')
% figure,plot(x,Dl,'-og')
% hold on
% plot(x,Du,'-or')
% xlabel('Number of defectives observed')
% ylabel('Lower and upper confidence bound of number of defectives in population')
% title('Two-sided symmetric confidence bounds as a function of number of defectives observed (N=1000, n=50)')
% hold off
% figure,plot(D(1:end),p(1:end),'-ok')
% xlabel('Number of defectives in population')
% ylabel('Coverage probability')
% title('Coverage probabilities of symmetric two-sided confidence interval as a function of number of defectives observed (N=1000, n=50)')

%% method 4 - clopper-pearson interval (exact computation using inverse of beta cdf)
Du4 = floor(N*betainv(1-alpha1,x+1,n-x));
Dl4 = ceil(N*betainv(alpha2,x,n-x+1));
[x' Dl4' Du4' Du4'-Dl4']