%figure(2)
tic

targetsignal = csvread("volterra.csv");
targetsignal10 = targetsignal(1:end-1);
%targetsignal2 = csvread("targetsignal2.csv");
%targetsignal3 = csvread("targetsignal3.csv");

%inputsignal =csvread("inputsignalcheck.csv");
%inputsignal2 =csvread("Data/inputsignal.csv");plot
%test2 = test1(:,2);
%learningweights = csvread("learningweights.csv");
%test4 = csvread("outputsignal.csv");
%test5 = csvread("inputsignal.csv");
learningmatrix = csvread("learningmatrix.csv");
learningmatrix2 = csvread("learningmatrix2.csv");
learningmatrix3 = csvread("learningmatrix3.csv");

%learningmatrix = learningmatrix(:, 1:end-1);
%learningmatrix2 = learningmatrix(10000:410000, 1:end-1);
%learningmatrix3 = learningmatrix(410000:425000, 1:end-1);

%learningmatrix2(:, end+1) = ones(1, 200000);
%learningmatrix3(:, end+1) = ones(1,15000);

%targetsignal2 = targetsignal(10000:410000);
%targetsignal3 = targetsignal(410000:425000);
%learningweights2 = learningmatrix2\targetsignal2;

%OutputVector = legarningmatrix3*learningweights;
%OutputVector2 = learningmatrix3*learningweights2;
%plot(OutputVector2)
%hold on
%plot(targetsignal3)

%plot(OutputVector)
%hold on
%plot(OutputVector2)
%learningweightstest = learningmatrix\targetsignal10;
%OutputVector10 = learningmatrix*learningweightstest;
%figure(10);
%plot(OutputVector10)




%learningmatrix*(learningmatrix\targetsignal)

%[m,n] = size(learningmatrix);
% = 1 + ((learninpgmatrix(:,1)-mean(learningmatrix(:,1))*(std(learningmatrix(:,1)/1))))x;

%newmatrix = (learningmatrix3 - repmat(mean(learningmatrix3), size(learningmatrix3,1), 1)) ./ repmat(std(learningmatrix3), size(learningmatrix3,1), 1);
newmatrix2 = (learningmatrix - repmat(mean(learningmatrix), size(learningmatrix,1), 1)) ./ repmat(std(learningmatrix), size(learningmatrix,1), 1);
newmatrix3 = (learningmatrix - repmat(mean(learningmatrix), size(learningmatrix,1), 1));
%plot(OutputMatrix, "Color", "Blue")
%hold on
%plot(targetsignal)

%chaoscheck = csvread("chaoscheck.csv");
%chaoscheck1 = csvread("chaoscheck1.csv");
%chaoscheck2 = csvread("chaoscheck2.csv");
%chaoscheck3 = csvread("chaoscheck3.csv");
%chaoscheck4 = csvread("chaoscheck4.csv");
%chaoscheck5 = csvread("chaoscheck5.csv");
%chaoscheck6 = csvread("chaoscheck6.csv");
%chaoscheck9 = csvread("chaoscheck9.csv");
%chaoscheck10 = csvread("chaoscheck10.csv");


%newmatrix = (learningmatrix - repmat(mean(learningmatrix), size(learningmatrix,1), 1)) ./ repmat(std(learningmatrix), size(learningmatrix,1), 1);
%plot(OutputMatrix, "Color", "Blue")
%hold on
%plot(targetsignal)

%plot(chaoscheck)
%hold on 
%plot(chaoscheck1)
%hold on
%plot(chaoscheck2)
%hold on
%plot(chaoscheck3);
%hold on 
%plot(chaoscheck4);
%hold on
%plot(chaoscheck5);
%
%plot(chaoscheck6);
%%plot(chaoscheck9);
%hold on
%plot(chaoscheck10);
%hold on 
%title("input force = 1.000 to input force  1.01")
%xlabel("Timesteps");
%ylabel("Sum of L matrix for all springs");

%NodePositionsx = csvread("NodePositionsx.csv");
%NodeVelocitiesx = csvread("NodeVelocitiesx.csv");
%NodeAccelerationsx = csvread("NodeAccelerationsx.csv");

%NodePositionsy = csvread("NodePositionsy.csv");
%NodeVelocitiesy = csvread("NodeVelocitiesy.csv");
%NodeAccelerationsy = csvread("NodeAccelerationsy.csv");


%plot(NodeAccelerationsx)
%hold on
%plot(NodeAccelerationsy)





%figure(3)
%stem(learningweights);
%hold on
%plot(test1)

%figure(1)
%plot((test5))
%learning_Volterra
%plot(normc(ImpulseResponse))
toc



%plot(test3)
%stem(test3);

%plot(test2(:,1),test2(:,2),test4(:,1),test4(:,2))


