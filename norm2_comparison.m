%------------------------------------------------------------------------------
% Filename: norm2_comparison.m
% 
% To determine errors in approximation of hypergeometric distribution
%
% This code takes population size and sample size to generate sampling
% distributions. Approximation errors are determined for normal and binomial
% approximations.
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
D = 0:N; %defectives in population
n = 20; %sample size

p = D/N; %proportion of defectives in population
x = 0:n; %defectives in sample
errnorm2 = zeros(N+1,2);
errbino = zeros(N+1,2);

for i = D
    zhype = hygepdf(x,N,i,n);
    znorm2 = normcdf(x+0.5,n*p(i+1),sqrt(n*p(i+1)*(1-p(i+1))*(N-n)/(N-1)))-normcdf(x-0.5,n*p(i+1),sqrt(n*p(i+1)*(1-p(i+1))*(N-n)/(N-1)));
    zbino = binopdf(x,n,p(i+1));

    yhype = hygecdf(x,N,D(i+1),n);
    ynorm2 = normcdf(x+0.5,n*p(i+1),sqrt(n*p(i+1)*(1-p(i+1))*(N-n)/(N-1)));
    ybino = binocdf(x,n,p(i+1));

    errnorm2(i+1,:) = [norm(zhype-znorm2,2) norm(yhype-ynorm2,2)];
    errbino(i+1,:) = [norm(zhype-zbino,2) norm(yhype-ybino,2)];
    
end

% plot(D,errnorm2(:,1),'-xr');
% hold on
% plot(D,errpois(:,1),'-+b');
% xlabel('Number of defectives in population');
% ylabel('Error');
% legend('normal (hyper)','poisson');
% str = sprintf('L2 norm of error in approximation of pdf for N = %d, n = %d',N,n);
% title(str);
% hold off

figure,plot(D,errnorm2(:,2),'-xr');
hold on
plot(D,errbino(:,2),'-+b');
xlabel('Number of defectives in population');
ylabel('Error');
legend('normal (hyper)','binomial');
str = sprintf('L2 norm of error in approximation of cdf for N = %d, n = %d',N,n);
title(str);
hold off