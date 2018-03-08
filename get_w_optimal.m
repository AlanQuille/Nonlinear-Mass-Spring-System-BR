
clear all
close all

L = csvread('learningmatrix.csv');
o = csvread('volterra.csv');

i = csvread('input.csv');
% o = i.*2;




washout = 100000;


L_ = L(washout+1:end,:);
% L_ones = ones(lenght(L_,1);

len = length(L_);
%L_full = [L;ones(len)];
o_ = o(washout+1:washout+len,:);
w  = L_\o_;
figure;stem(w)

y = L_*w;
figure;plot(y,'b');
hold on;plot(o_,'r--');
% plotting
% L_= mapstd(L_');
% figure;plot(L_');
mean((y-o_).^2)


L2 = mapstd(L');
figure;plot(L2');


