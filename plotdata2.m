%figure(2)
tic
test1 = csvread("targetsignal.csv");
inputsignal =csvread("inputsignalcheck.csv");
test2 = test1(:,2);
test3 = csvread("learningweights.csv");
%test4 = csvread("outputsignal.csv");
%test5 = csvread("inputsignal.csv");
test5 = csvread("learningmatrix.csv");
test6 = pinv(test5);
test3 = test6 * test2;
test7 = test5 * test3;

%plot(test7);
%hold on
%plot(test1)

figure(1)
plot((test5))
%learning_Volterra
%plot(normc(ImpulseResponse))
toc



%plot(test3)
%stem(test3);

%plot(test2(:,1),test2(:,2),test4(:,1),test4(:,2))


