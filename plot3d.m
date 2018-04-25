matrix = csvread("k_and_d_and_MSE.csv");
matrix2 = matrix(:, 1:2);
matrix3 = matrix(:, 4);
matrix5 = [matrix2 matrix3];

x1 = matrix5(:,1:1);
y1 = matrix5(:,2:2);
z1 = matrix5(:,3:3);
plot3(x1, y1, z1)