figure(2)
test2 = csvread("targetsignal.csv");
test3 = csvread("learningweights.csv");
test4 = csvread("outputsignal.csv");
test5 = csvread("inputsignal.csv");

%plot(test5)
%stem(test3);

plot(test2(:,1),test2(:,2),test4(:,1),test4(:,2))


