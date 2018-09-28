MSE=load("MSE_and_scaling_factor_1.csv")
MSE2 = MSE;
minMSE = min(MSE);
MSE = MSE(:,2)

%plot(MSE(:,1), MSE(:,2))
%xticklabels(MSE(:,1))
%MSE=(MSE-MSE)./MSE;
plot(MSE)
xticks(0:20)
%OneD10DegDiff=(OneD10Deg-OneDBaseline)./OneDBaseline;
xticklabels([0 10^(-19) 10^(-18) 10^(-17) 10^(-16) 10^(-15) 10^(-14) 10^(-13) 10^(-12) 10^(-11) 10^(-10) 10^(-9) 10^(-8) 10^(-7) 10^(-6) 10^(-5) 10^(-4) 10^(-3) 10^(-2) 10^(-1) 1])
%ylim([0 1000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000])
%ylim([0,10]);
xlabel({'Scaling Factor T','s'})
ylabel({'MSE'})