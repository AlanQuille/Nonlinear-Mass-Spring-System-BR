figure(1);
axis([0 10 0 10])
for i=1:999
s = num2str(i*0.001);
if(length(s)==4)
    s =strcat(s, "0");
end

if(length(s)==3)
    s = strcat(s, "00");
end
s = strcat(s, ".csv");
disp(s)
test = csvread(s);
x = [test(:,1),test(:,3)];
y = [test(:,2),test(:,4)];
z =  sin(i*0.001);
plot(x,y,'-o','MarkerSize',10, 'MarkerEdgeColor','red');
drawnow;
end