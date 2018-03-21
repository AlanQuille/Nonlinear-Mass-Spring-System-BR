% read and plot L matrix
close all
clear all

addpath('');

L = csvread('learningmatrix.csv');
% w = csvread('outputweights.csv');
o = csvread('outputsignal.csv');

L_ = mapstd(L');
figure;plot(L_');
figure;stem(w);

figure;plot(o(:,2))