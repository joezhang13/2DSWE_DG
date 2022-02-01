clear;

%% Free-stream and free-stream preservation test
mesh = readgri('city0.gri');
CFL0 = 0.95;
normR1 = zeros(1, 3);
dt = zeros(1, 3);
for p = 0 : 2
    Nstep = 1;
    [normR1(p + 1), ~] = FreeStream(mesh, CFL0, Nstep, p);
    Nstep = 1500;
    [normR2, dt(p + 1)] = FreeStream(mesh, CFL0, Nstep, p);
    figure;
    plot(normR2, 'LineWidth', 1.25);
    xlabel('Number of steps');
    ylabel('L1 norm of residual');
    set(gca, 'FontSize', 12);
end
disp(['Free-stream test (p = ', num2str(0:2), '):']);
disp(['L1 norm is ', num2str(normR1)]);
disp(['Time step is ', num2str(dt)]);

%% Forces on the buildings
CFL0 = 0.95; T = 2;
mesh0 = readgri('city0.gri');
mesh1 = readgri('city1.gri');
Force0 = cell(4, 1); t0 = cell(4, 1);
Force1 = cell(4, 1); t1 = cell(4, 1);
for p = 0 : 3
    [U0, dt0] = initial(mesh0, CFL0, p);
    [~, Force0{p+1}] = DGSolver(mesh0, U0, dt0, T, p);
    t0{p+1} = linspace(0, T, ceil(T/dt0) + 1);
    [U1, dt1] = initial(mesh1, CFL0, p);
    [~, Force1{p+1}] = DGSolver(mesh1, U1, dt1, T, p);
    t1{p+1} = linspace(0, T, ceil(T/dt1) + 1);
end

% load forces.mat;    %read the simulation results

for i = 1 : 4
    figure;
    subplot(1, 2, 1);
    for p = 0 : 3
        plot(t0{p+1}, Force0{p+1}(:, 1, i),'LineWidth', 1.5);
        hold on
    end
    xlabel('t');
    ylabel('F_x');
    legend('p = 0', 'p = 1', 'p = 2', 'p = 3');
    set(gca, 'FontSize', 12);
    subplot(1, 2, 2);
    for p = 0 : 3
        plot(t0{p+1}, Force0{p+1}(:, 2, i),'LineWidth', 1.5);
        hold on
    end
    xlabel('t');
    ylabel('F_y');
    legend('p = 0', 'p = 1', 'p = 2', 'p = 3');
    set(gca, 'FontSize', 12);
%     sgtitle(['Forces on building ', num2str(i), ' for the coarse mesh']);
end
for i = 1 : 4
    figure;
    subplot(1, 2, 1);
    for p = 0 : 3
        plot(t1{p+1}, Force1{p+1}(:, 1, i),'LineWidth', 1.5);
        hold on
    end
    xlabel('t');
    ylabel('F_x');
    legend('p = 0', 'p = 1', 'p = 2', 'p = 3');
    set(gca, 'FontSize', 12);
    subplot(1, 2, 2);
    for p = 0 : 3
        plot(t1{p+1}, Force1{p+1}(:, 2, i),'LineWidth', 1.5);
        hold on
    end
    xlabel('t');
    ylabel('F_y');
    legend('p = 0', 'p = 1', 'p = 2', 'p = 3');
    set(gca, 'FontSize', 12);
%     sgtitle(['Forces on building ', num2str(i), ' for the fine mesh']);
end

%% Contour plot of water height
CFL0 = 0.95;
mesh0 = readgri('city0.gri');
T = [0, 0.05, 0.1, 0.15, 0.2, 0.25];
U = cell(4, 1);
for p = 0 : 3
    [U0, dt0] = initial(mesh0, CFL0, p);
    [U{p+1}, ~] = DGSolver(mesh0, U0, dt0, T, p);
end

% load solutions.mat;                  %read the simulation results

for i = 1 : 6
    for p = 0 : 3
        solplot(mesh0, U{p+1}(:, :, i), p, T(i));
    end
end
