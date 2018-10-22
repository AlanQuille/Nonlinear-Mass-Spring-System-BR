stable = load("buffer_zone_stable_3.csv");
unstable = load("buffer_zone_unstable_6_m01.csv");
figure(1);
scatter(stable(:,1), stable(:,2), 50, 'g');
hold on
scatter(unstable(:,1), unstable(:,2), 50, 'r', 'd');
xlabel('k_{1}')
ylabel('k_{3}')
legend('stable','unstable')
figure(2);
scatter(stable(:,1), stable(:,3), 50, 'g');
hold on
scatter(unstable(:,1), unstable(:,3), 50, 'r', 'd');
xlabel('k_{1}')
ylabel('d_{1}')
legend('stable','unstable')
figure(3);
scatter(stable(:,1), stable(:,4), 50, 'g');
hold on
scatter(unstable(:,1), unstable(:,4), 50, 'r', 'd');
xlabel('k_{1}')
ylabel('d_{3}')
legend('stable','unstable')
figure(4);
scatter(stable(:,2), stable(:,3), 50, 'g');
hold on
scatter(unstable(:,2), unstable(:,3), 50, 'r', 'd');
xlabel('k_{3}')
ylabel('d_{1}')
legend('stable','unstable')
figure(5);
scatter(stable(:,2), stable(:,4), 50, 'g');
hold on
scatter(unstable(:,2), unstable(:,4), 50, 'r', 'd');
xlabel('k_{3}')
ylabel('d_{3}')
legend('stable','unstable')
figure(6);
scatter(stable(:,3), stable(:,4), 50, 'g');
hold on
scatter(unstable(:,3), unstable(:,4), 50, 'r', 'd');
xlabel('d_{1}')
ylabel('d_{3}')
legend('stable','unstable')
