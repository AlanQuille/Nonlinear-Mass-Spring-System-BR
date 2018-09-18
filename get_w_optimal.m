
clear all
close all

addpath('');

%L = csvread('learningmatrix.csv');
%L = L(:, 1:end-1);
%L = csvread('learningmatrix3.csv');
%L = L(:, 1:end-1);
inputsig = csvread('Data/inputsignal.csv');
targetsig = csvread('Data/narma.csv');
outputsignal = csvread('1_outputsignal.csv');
%o = csvread('targetsignal.csv');

figure(2)
plot(targetsig(403000:418000));
hold on
plot(outputsignal);

%i = csvread('input.csv');
% o = i.*2;




%washout = 20000;


%L_ = L(washout+200000:end,:);
%L_ = L;
%L_ones = ones(length(L_,1));

%len = length(L_);
%L_full = [L;ones(len)];
%o_ = o(washout+1:washout+len,:);
%o_ = o;y
%w1  = L_\o_;
%TRY TRY TRY TRY TRY TRY TRY TRY
%w = csvread("learningweights.csv");
%figure;stem(w)

%y = L_*w;
%figure;plot(y,'b');
%hold on;plot(o_,'r--');
% plotting
% L_= mapstd(L_');
% figure;plot(L_');
%mean((y-o_).^2)


%L2 = mapstd(L');
%figure;plot(L2');


