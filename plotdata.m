
for i=1:1
s = num2str(i*0.001);
if(length(s)==4)
    s =strcat(s, "0");
end

if(length(s)==3)
    s = strcat(s, "00");
end
s1 = strcat(s, "X.csv");
s2 = strcat(s, "Y.csv");
s3 = strcat(s, "Z.csv");
%disp(s)

test = csvread(s1)
test2 = csvread(s2);
test3 = csvread(s3);
x = test;
y = test2;
z = test3;


%x = [5,3,1];c
%y = [1,4,2];
%plot(x,y,'-o','MarkerSize',10, 'MarkerEdgeColor','red');
%plot(x,y,'Layout','force')


%x=[x_pos(nodea), x_pos(nodeb), x_pos(nodea), x_pos(nodeb)]
%y=[y_pos(nodea), y_pos(nodeb), y_pos(nodea)]

s = csvread("s.csv")
t = csvread("t.csv")
G = graph(s,t);

%G=addnode(G,1);









%x = [test(1,1) test(1,3) test(2,1)];
%y = [test(1,2) test(1,4) test(2,2)];

%test(:,1)
plot(G,'XData',x,'YData',y, 'ZData', z);
%xlim([-0.215 -0.19])
%ylim([0.4 0.45])
axis([-10 10 -10 10]);
pause(0.5);
drawnow;
end