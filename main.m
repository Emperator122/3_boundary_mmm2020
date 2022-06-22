clear;
addpath("Functions");

alpha = 0.75;
sigma = 1-alpha./2;

% точное решение
u_exact = @(hi,t)(t.^3+3.*t.^2+1).*(hi.^6+2*hi+3); 

% Начальное условие
initialCondition = @(hi)hi.^6+2*hi+3;

% Границы для t
a_T = 0;
b_T = 1;
M = 40; % Колво разбиений отрезка [a_T; b_T]

% Границы для hi
a_x = 0;
b_x = 1;
N = 40; % количества разбиений каждого отрезка

% Шаги
tau = (b_T-a_T)/M;
h = (b_x-a_x)/N;

t = a_T:tau:b_T;
hi = a_x:h:b_x;

% Все про i_0
i_0 = N-2; % выбор самой точки
x_i_0 = hi(i_0); 
x_i_0_plus_1 = hi(i_0+1);
x_0 = x_i_0+h./2; % задание узла
x_left_i_0 = (x_i_0_plus_1-x_0)./h;
x_right_i_0 = (x_0-x_i_0)./h;

% для третьего краевого
beta_1 = @(t)3+t.^2;
beta_2 = @(t)2.*(2-sin(t));
mu_1 = @(t)(t.^3+3.*t.^2+1).*(5+3.*t.^2);
mu_2 = @(t)(t.^3+3.*t.^2+1).*(2-sin(t)).*20;

% k(hi,t)
k = @(hi, t)2-sin(hi.*t);

%q(hi,t)
q = @(hi,t)1-cos(hi.*t);

% f(hi,t)
f = @(hi,t) ( 6.*t.^(3-alpha)./gamma(4-alpha) + ...
   6.*t.^(2-alpha)./gamma(3-alpha)).*(hi.^6+2.*hi+3)...
   -(t.^3+3.*t^2+1) .* ( -t.*cos(hi.*t).*(6.*hi.^5+2) ...
   + (2-sin(hi.*t)).*30.*hi.^4 ...
   - (x_0.^6+2.*x_0+3).*(1-cos(hi.*t)) );

% параметры графика
layer = length(t);

% Расчеты
tic
u = zeros(N+1, M+1); % первая для hi, вторая для t
u(:, 1) = initialCondition(hi);

for n = 1:M
    % зададим матрицу и столбец свободных членов для метода прогонки
    matrix = zeros(N+1, N+1);
    matrix_f = zeros(N+1, 1);

    % необходимые значения
    step_coef = sigma./h.^2;
    sigmaStep = t(n)+sigma.*tau;
    gamma_value = tau.^(-alpha) ./ gamma(2-alpha);

    % вычислим "c"
    c_values = calculate_c(n, alpha, sigma);

    unknown_coef = kaputo_der_unknown_coef(c_values, gamma_value);
    % загоним в матрицу краевые условия
    mu_1_wave = mu_1(sigmaStep) + 0.5.*h.*f(hi(1), sigmaStep);
    known_part = kaputo_der_known_part(c_values, 1, n, u, gamma_value);
    matrix(1, 1) = -0.5.*h.*unknown_coef ...
        - sigma.*beta_1(sigmaStep) - sigma./h.*k(hi(2)-h./2, sigmaStep);
    matrix(1, 2) = sigma./h.*k(hi(2)-h./2, sigmaStep);
    matrix(1, i_0) = matrix(1, i_0) + ...
        -0.5.*h.*sigma.*q(hi(1), sigmaStep).*x_left_i_0;
    matrix(1, i_0+1) = matrix(1, i_0+1) + ...
        -0.5.*h.*sigma.*q(hi(1), sigmaStep).*x_right_i_0;

    matrix_f(1) = 0.5.*h.*known_part ...
        + (1-sigma).*beta_1(sigmaStep).*u(1, n) - mu_1_wave ...
        + 0.5.*h.*(1-sigma).*q(hi(1), sigmaStep).*...
        (x_left_i_0.*u(i_0, n) + x_right_i_0.*u(i_0+1, n)) - ...
        (1-sigma).*(u(2, n)-u(1, n))./h.*k(hi(2)-h./2, sigmaStep);
% 
    mu_2_wave = mu_2(sigmaStep) + 0.5.*h.*f(hi(N+1), sigmaStep);
    known_part = kaputo_der_known_part(c_values, N+1, n, u, gamma_value);
    matrix(N+1,N) = sigma./h.*k(hi(N+1)-h./2, sigmaStep);
    matrix(N+1,N+1) = -0.5.*h.*unknown_coef ...
        - sigma.*beta_2(sigmaStep) - sigma./h.*k(hi(N+1)-h./2, sigmaStep);
    matrix(N+1, i_0) = matrix(N+1, i_0) + ...
        -0.5.*h.*sigma.*q(hi(N+1), sigmaStep).*x_left_i_0;
    matrix(N+1, i_0+1) = matrix(N+1, i_0+1) + ...
        -0.5.*h.*sigma.*q(hi(N+1), sigmaStep).*x_right_i_0;

    matrix_f(N+1) = 0.5.*h.*known_part ...
        + (1-sigma).*beta_2(sigmaStep).*u(N+1, n) - mu_2_wave ...
        + 0.5.*h.*(1-sigma).*q(hi(N+1), sigmaStep).*...
        (x_left_i_0.*u(i_0, n) + x_right_i_0.*u(i_0+1, n)) + ...
        (1-sigma).*(u(N+1, n)-u(N, n))./h.*k(hi(N+1)-h./2, sigmaStep);

    % загоним в матрицу остальные точки
    for i = 2:N
        known_part = kaputo_der_known_part(c_values, i, n, u, gamma_value);

        % для начала заполним саму матрицу
        matrix(i, i-1) = step_coef.*...
            k(hi(i)-h./2, sigmaStep);
        matrix(i, i) = -1.*(unknown_coef + ...
            (k(hi(i+1)-h./2, sigmaStep)+k(hi(i)-h./2, sigmaStep)).*step_coef);
        matrix(i, i+1) = step_coef.*k(hi(i+1)-h./2, sigmaStep);
        matrix(i, i_0) = matrix(i, i_0) + ...
            -1.*sigma.*q(hi(i), sigmaStep).*x_left_i_0;
        matrix(i, i_0+1) = matrix(i, i_0+1) + ... 
            -1.*sigma.*q(hi(i), sigmaStep).*x_right_i_0;
        
        % а затем столбец свободных членов
        % для заполнения столбца с.ч. необходимо вычислять сумму
        % вычислим ее
        sum = known_part;
        sum = sum - f(hi(i), sigmaStep);
        sum = sum - (1-sigma).*(k(hi(i+1)-h./2, sigmaStep).*u(i+1, n) ...
            -(k(hi(i+1)-h./2, sigmaStep)+k(hi(i)-h./2, sigmaStep)).*u(i, n)...
            +k(hi(i)-h./2, sigmaStep).*u(i-1, n))./h.^2;
        sum = sum + (1-sigma).*u(i_0, n).*q(hi(i), sigmaStep).*x_left_i_0;
        sum = sum + (1-sigma).*u(i_0+1, n).*q(hi(i), sigmaStep).*x_right_i_0;
        
        % результат вычислений закидываем в столбец с.ч.
        matrix_f(i) = sum;
    end
    % находим корни и записываем их на новый слой
    roots = matrix\matrix_f;
    u(:, n+1) = roots;
end
toc
delta = abs(u_exact(hi', t)-u);
plot(hi, u(:, layer), '-r+', hi, u_exact(hi', t(layer)), '-bo');
legend("Численное решение", "Точное решение");
max(max(delta))


function known_part = kaputo_der_known_part(c_values, i, j, u, gamma_value)
        known_part = (sum(c_values(2:j) .* (u(i, j:-1:2)-u(i, j-1:-1:1)))-c_values(1).*u(i, j)) ...
        .* gamma_value;    
end

function unknown_coef = kaputo_der_unknown_coef(c_values, gamma_value)
    unknown_coef = c_values(1).* gamma_value;
end