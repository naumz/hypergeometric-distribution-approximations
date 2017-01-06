%------------------------------------------------------------------------------
% Filename: approx_comparison.m
% 
% To determine probability distributions of various approximations to
% hypergeometric distribution.
%		
% This code takes population size, number of defectives and sample size to
% generate sampling distributions. The hypergeometric distribution (exact) is
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

N = 1000; %population size
D = 250; %defectives in population
n = 20; %sample size

p = D/N; %proportion of defectives in population
x = 0:n; %defectives in sample

zhype = hygepdf(x,N,D,n);
zbino = binopdf(x,n,p);
znorm = normcdf(x+0.5,n*p,sqrt(n*p*(1-p)))-normcdf(x-0.5,n*p,sqrt(n*p*(1-p)));
znorm2 = normcdf(x+0.5,n*p,sqrt(n*p*(1-p)*(N-n)/(N-1)))-normcdf(x-0.5,n*p,sqrt(n*p*(1-p)*(N-n)/(N-1)));
zpois = poisspdf(x,n*p);

yhype = hygecdf(x,N,D,n);
ybino = binocdf(x,n,p);
ynorm = normcdf(x+0.5,n*p,sqrt(n*p*(1-p)));
ynorm2 = normcdf(x+0.5,n*p,sqrt(n*p*(1-p)*(N-n)/(N-1)));
ypois = poisscdf(x,n*p);

errbino = [norm(yhype-ybino,2) norm(zhype-zbino,2)]
errnorm = [norm(yhype-ynorm,2) norm(zhype-znorm,2)]
errnorm2 = [norm(yhype-ynorm2,2) norm(zhype-znorm2,2)]
errpois = [norm(yhype-ypois,2) norm(zhype-zpois,2)]

plot(x,zhype,'-+g');
hold on
plot(x,zbino,'-ob');
plot(x,znorm,'-xr');
plot(x,znorm2,'-dm');
plot(x,zpois,'-*k');
xlabel('Number of defectives observed in sample');
ylabel('pdf')
str = sprintf('pdf of different distribtutions for N = %d, D = %d, n = %d',N,D,n);
title(str);
legend('hypergeometric','binomial','normal (bino)','normal (hyper)','poisson');
hold off

figure,plot(x,yhype,'-+g');
hold on
plot(x,ybino,'-ob');
plot(x,ynorm,'-xr');
plot(x,ynorm2,'-dm');
plot(x,ypois,'-*k');
xlabel('Number of defectives observed in sample');
ylabel('cdf')
str = sprintf('cdf of different distribtutions for N = %d, D = %d, n = %d',N,D,n);
title(str);
legend('hypergeometric','binomial','normal (bino)','normal (hyper)','poisson');
hold off