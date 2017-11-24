function plotdata()

%maxtimesteps = 10/0.001;
xlim([0,10]);
ylim([0,10]);
%string nodes.csv;
%for i = 0:maxtimesteps-1
    test = csvread("nodes.csv");
    %insertBefore(string,".",i)
    %x=[test(:,1+i) test(:,3+i)];
    %y=[test(:,2+i) test(:,4+i)];
    x=[test(:,1) test(:,3)];
    y=[test(:,2) test(:,4)];
    disp(x);
    plot(x,y,'Marker', 'o')
end 
