%figure(2)
tic

targetsignal = csvread("targetsignal.csv");
inputsignal =csvread("inputsignalcheck.csv");
%inputsignal2 =csvread("Data/inputsignal.csv");
%test2 = test1(:,2);
%learningweights = csvread("learningweights.csv");
%test4 = csvread("outputsignal.csv");
%test5 = csvread("inputsignal.csv");
learningmatrix = csvread("learningmatrix.csv");
learningweights = learningmatrix\targetsignal;
OutputMatrix = learningmatrix*learningweights;

%[m,n] = size(learningmatrix);
% = 1 + ((learningmatrix(:,1)-mean(learningmatrix(:,1))*(std(learningmatrix(:,1)/1))))x;

newmatrix = (learningmatrix - repmat(mean(learningmatrix), size(learningmatrix,1), 1)) ./ repmat(std(learningmatrix), size(learningmatrix,1), 1);
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

NodePositionsx = csvread("NodePositionsx.csv");
NodeVelocitiesx = csvread("NodeVelocitiesx.csv");
NodeAccelerationsx = csvread("NodeAccelerationsx.csv");

NodePositionsy = csvread("NodePositionsy.csv");
NodeVelocitiesy = csvread("NodeVelocitiesy.csv");
NodeAccelerationsy = csvread("NodeAccelerationsy.csv");


plot(NodeAccelerationsx)
hold on
plot(NodeAccelerationsy)





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


