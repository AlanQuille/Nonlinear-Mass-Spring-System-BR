figure(2)
test2 = csvread("targetsignal.csv");
test3 = csvread("learningweights.csv");
test4 = csvread("outputsignal.csv");
test5 = csvread("newinput.csv");
test6 = csvread("impulseresponse.csv");
test6 = normc(test6);
length(test6)
stem(test3);

%plot(test5)
%plot(test6)

%plot(test6)
%plot(test4(:,1), test4(:,2), test2(:,1), test2(:,2))
plot(test2(:,1), test2(:,2))
%plot(test4(:,1), test4(:,2))


mean(test6)
