#include <iostream>
#include <iomanip>
#include <cmath>
#include <limits>
#include <list>

// Test 
double test_system(double x, double t)
{
    //return -50 * (x -cos(t));
    return  -x * x * x +5 * sin(10 * t);
}


// Options for Simulation environment
struct simulation_options
{
    const double t_start = 0.0;
    const double t_end = 1.0;
    const double h = 0.05;
};

// Options for implicit solver
struct solver_options
{
    const int max_iter = 2;
    std::list<int> iter_result;

    const double rel_tol = 1e-8;
    std::list<double> rel_tol_result;

    const double abs_tol = 1e-6;
    std::list<double> abs_tol_result;

    const int order = 3;
};


// Newton iteration 
double newton(double x, double t, double h, double x_i, double (*test_system)(double, double), solver_options& options, double dfdx, bool reset_counter = true)
{
    //static const double eps = std::pow(std::numeric_limits<double>::epsilon(), 1.0/2.0);

    //double dfdx = (test_system(x_i + eps / 2.0, t+h) - test_system(x_i - eps / 2.0, t+h)) / eps;
    //double dfdx = (test_system(x + eps / 2.0, t+h) - test_system(x - eps / 2.0, t+h)) / eps;
    double dg = 1 - h * dfdx;

    double g = x_i - x - h * test_system(x_i, t+h);

    double x_i_p_1 = x_i - g / dg;

    static int counter = 1;
    if (reset_counter)
        counter = 1;
    else
        ++counter;

    double abs_tol = std::abs(x_i - x_i_p_1);
    double rel_tol = std::abs(1 - x_i / x_i_p_1);

    if ((abs_tol <= options.abs_tol) ||
        (rel_tol <= options.rel_tol) ||
        (counter >= options.max_iter))
    {
        options.iter_result.push_back(counter);
        options.abs_tol_result.push_back(abs_tol);
        options.rel_tol_result.push_back(rel_tol);

        return x_i_p_1;
    }
    else
        return newton(x, t, h, x_i_p_1, test_system, options, dfdx, false);
}

// Euler
double ode1(double x, double t, double h, double (*test_system)(double, double))
{
    double k1 = test_system(x, t);
    return x + h * k1;
}

// Heun
double ode2(double x, double t, double h, double (*test_system)(double, double))
{
    double k1 = test_system(x, t);
    double k2 = test_system(x + h * k1, t + h);
    return x + 1.0/2.0 * h * (k1 + k2);
}

// Bogacki–Shampine
double ode3(double x, double t, double h, double (*test_system)(double, double))
{
    double k1 = test_system(x, t);
    double k2 = test_system(x + 1.0/2.0 * h * k1, t + 1.0/2.0 * h);
    double k3 = test_system(x + 3.0/4.0 * h * k2, t + 3.0/4.0 * h);
    return x + h * (2.0/9.0 * k1 + 1.0/3.0 * k2 + 4.0/9.0 * k3);
}

// Runge-Kutta 4
double ode4(double x, double t, double h, double (*test_system)(double, double))
{
    double k1 = test_system(x, t);
    double k2 = test_system(x + 1.0/2.0 * h * k1, t + 1.0/2.0 * h);
    double k3 = test_system(x + 1.0/2.0 * h * k2, t + 1.0/2.0 * h);
    double k4 = test_system(x + h * k3, t + h);
    return x + h * 1.0/6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
}

// Dormand-Prince (fifth order)
double ode5(double x, double t, double h, double (*test_system)(double, double))
{
    double k1 = test_system(x, t);
    double k2 = test_system(x + 1.0/5.0 * h * k1, t + 1.0/5.0 * h);
    double k3 = test_system(x + 3.0/40.0 * h * k1 + 9.0/40.0 * h * k2, t + 3.0/10.0 * h);
    double k4 = test_system(x + 44.0/45.0 * h * k1 - 56.0/15.0 * h * k2 + 32.0/9.0 * h * k3, t + 8.0/10.0 * h);
    double k5 = test_system(x + 19372.0/6561.0 * h * k1 - 25360.0/2187.0 * h * k2 + 64448.0/6561.0 * h * k3 - 212.0/729.0 * h * k4, t + 8.0/9.0 * h);
    double k6 = test_system(x + 9017.0/3168.0 * h * k1 - 355.0/33.0 * h * k2 + 46732.0/5247.0 * h * k3 + 49.0/176.0 * h * k4 - 5103.0/18656.0 * h * k5, t + h);
    return x + h * (35.0/384.0 * k1 + 500.0/1113.0 * k3 + 125.0/192.0 * k4 - 2187.0/6784.0 * k5 + 11.0/84.0 * k6);
}

// Dormand-Prince (eight order)
double ode8(double x, double t, double h, double (*test_system)(double, double))
{
    double k1 = test_system(x, t);
    double k2 = test_system(x + 1.0/18.0 * h * k1, t + 1.0/18.0 * h);
    double k3 = test_system(x + 1.0/48.0 * h * k1 + 1.0/16.0 * h * k2, t + 1.0/12.0 * h);
    double k4 = test_system(x + 1.0/32.0 * h * k1 + 3.0/32.0  * h * k3, t + 1.0/8.0 * h);
    double k5 = test_system(x + 5.0/16.0 * h * k1 - 75.0/64.0 * h * k3 + 75.0/64.0 * h * k4, t + 5.0/16.0 * h);
    double k6 = test_system(x + 3.0/80.0 * h * k1 + 3.0/16.0 * h * k4 + 3.0/20.0 * h * k5, t + 3.0/8.0 * h);
    double k7 = test_system(x + 29443841.0/614563906.0 * h * k1 + 77736538.0/692538347.0 * h * k4 - 28693883.0/1125000000.0 * h * k5 + 23124283.0/1800000000.0 * h * k6, t + 59.0/400.0 * h);
    double k8 = test_system(x + 16016141.0/946692911.0 * h * k1 + 61564180.0/158732637.0 * h * k4 + 22789713.0/633445777.0 * h * k5 + 545815736.0/2771057229.0 * h * k6 - 180193667.0/1043307555.0 * h * k7, t + 93.0/200.0 * h);
    double k9 = test_system(x + 39632708.0/573591083.0 * h * k1 - 433636366.0/683701615.0 * h * k4 - 421739975.0/2616292301.0 * h * k5 + 100302831.0/723423059.0 * h * k6 + 790204164.0/839813087.0 * h * k7 + 800635310.0/3783071287.0 * h * k8, t + 5490023248.0/9719169821.0 * h);
    double k10 = test_system(x + 246121993.0/1340847787.0 * h * k1 - 37695042795.0/15268766246.0 * h * k4 - 309121744.0/1061227803.0 * h * k5 - 12992083.0/490766935.0  * h * k6 + 6005943493.0/2108947869.0  * h * k7 + 393006217.0/1396673457.0 * h * k8 + 123872331.0/1001029789.0 * h * k9, t + 13.0/20.0 * h);
    double k11 = test_system(x - 1028468189.0/846180014.0 * h * k1 + 8478235783.0/508512852.0 * h * k4 + 1311729495.0/1432422823.0 * h * k5 - 10304129995.0/1701304382.0 * h * k6 - 48777925059.0/3047939560.0 * h * k7 + 15336726248.0/1032824649.0 * h * k8 - 45442868181.0/3398467696.0 * h * k9 + 3065993473.0/597172653.0  * h * k10, t + 1201146811.0/1299019798.0 * h);
    double k12 = test_system(x + 185892177.0/718116043.0 * h * k1 - 3185094517.0/667107341.0 * h * k4 - 477755414.0/1098053517.0 * h * k5 - 703635378.0/230739211.0 * h * k6 + 5731566787.0/1027545527.0 * h * k7 + 5232866602.0/850066563.0 * h * k8 - 4093664535.0/808688257.0 * h * k9 + 3962137247.0/1805957418.0 * h * k10 + 65686358.0/487910083.0 * h * k11, t + h);
    double k13 = test_system(x + 403863854.0/491063109.0 * h * k1 - 5068492393.0/434740067.0 * h * k4 - 411421997.0/543043805.0 * h * k5 + 652783627.0/914296604.0 * h * k6 + 11173962825.0/925320556.0 * h * k7 - 13158990841.0/6184727034.0 * h * k8 + 3936647629.0/1978049680.0 * h * k9 - 160528059.0/685178525.0 * h * k10 + 248638103.0/1413531060.0 * h * k11, t + h);
    
    return x + h * (14005451.0/335480064.0 * k1 - 59238493.0/1068277825.0 * k6 + 181606767.0/758867731.0 * k7 + 561292985.0/797845732.0 * k8 - 1041891430.0/1371343529.0 * k9 + 760417239.0/1151165299.0 * k10 + 118820643.0/751138087.0 * k11 - 528747749.0/2220607170.0 * k12 + 1.0/4.0 * k13);
    //return x + h * (13451932.0/455176623.0 * k1 - 808719846.0/976000145.0 * k6 + 1757004468.0/5645159321.0 * k7 + 656045339.0/265891186.0 * k8 - 3867574721.0/1518517206.0 * k9 + 465885868.0/322736535.0 * k10 + 53011238.0/667516719.0 * k11 - 2.0/45.0 * k12);
}

// Backward Euler
double ode1be(double x, double t, double h, double (*test_system)(double, double), solver_options& options)
{
    static const double eps = std::pow(std::numeric_limits<double>::epsilon(), 1.0 / 2.0);
    double dfdx = (test_system(x + eps / 2.0, t) - test_system(x - eps / 2.0, t)) / eps;
    return newton(x, t, h, x, test_system, options, dfdx);
}

double ode14x(double x, double t, double h, double (*test_system)(double, double), solver_options& options)
{
    static const double eps = std::pow(std::numeric_limits<double>::epsilon(), 1.0 / 2.0);

    double dfdx = (test_system(x + eps / 2.0, t) - test_system(x - eps / 2.0, t)) / eps;

    auto newton_intermediate_steps = [](double x, double t, double h, double x_i, int number_of_steps, double (*test_system)(double, double), solver_options& options, double dfdx)
    {
        double hN = h / number_of_steps;
        double result = newton(x, t, hN, x, test_system, options, dfdx);

        for (size_t i = 1; i < number_of_steps; i++)
            result = newton(result, t + i * hN, hN, result, test_system, options, dfdx);

        return result;
    };

    double result1{ std::numeric_limits<double>::quiet_NaN() };
    double result2{ std::numeric_limits<double>::quiet_NaN() };
    double result3{ std::numeric_limits<double>::quiet_NaN() };
    double result4{ std::numeric_limits<double>::quiet_NaN() };

    if (options.order > 0)
        result1 = newton_intermediate_steps(x, t, h, x, 12, test_system, options, dfdx);

    if (options.order > 1)
        result2 = newton_intermediate_steps(x, t, h, x, 8, test_system, options, dfdx);

    if (options.order > 2)
        result3 = newton_intermediate_steps(x, t, h, x, 6, test_system, options, dfdx);

    if (options.order > 3)
        result4 = newton_intermediate_steps(x, t, h, x, 4, test_system, options, dfdx);

    if (options.order > 1)
        result1 = result2 - 3 * (result2 - result1);

    if (options.order > 2)
    {
        result2 = result3 - 4 * (result3 - result2);
        result1 = result2 - 2 * (result2 - result1);
    }

    if (options.order > 3)
    {
        result3 = result4 - 3 * (result4 - result3);
        result2 = result3 - 2 * (result3 - result2);
        result1 = result2 - 1.5 * (result2 - result1);
    }

    return result1;
}

int main()
{
    double x_ode1 = 4;
    double x_ode2 = x_ode1;
    double x_ode3 = x_ode1;
    double x_ode4 = x_ode1;
    double x_ode5 = x_ode1;
    double x_ode8 = x_ode1;
    double x_ode1i = x_ode1;
    double x_ode14x = x_ode1;

    simulation_options sim_options;

    std::cout << "t: " << std::setw(6) << sim_options.t_start << " x: " << std::setw(6) << x_ode1 << std::endl << std::endl;

    for (double t = sim_options.t_start; t < sim_options.t_end; t+= sim_options.h)
    {
        solver_options options_ode1i;
        solver_options options_ode14x;

        std::cout << "t: " << t+sim_options.h << std::endl;

        x_ode1 = ode1(x_ode1, t, sim_options.h, test_system);
        std::cout << "ode1      x: " << std::setw(6) << x_ode1 << std::endl;

        x_ode2 = ode2(x_ode2, t, sim_options.h, test_system);
        std::cout << "ode2      x: " << std::setw(6) << x_ode2 << std::endl;

        x_ode3 = ode3(x_ode3, t, sim_options.h, test_system);
        std::cout << "ode3      x: " << std::setw(6) << x_ode3 << std::endl;

        x_ode4 = ode4(x_ode4, t, sim_options.h, test_system);
        std::cout << "ode4      x: " << std::setw(6) << x_ode4 << std::endl;

        x_ode5 = ode5(x_ode5, t, sim_options.h, test_system);
        std::cout << "ode5      x: " << std::setw(6) << x_ode5 << std::endl;

        x_ode8 = ode8(x_ode8, t, sim_options.h, test_system);
        std::cout << "ode8      x: " << std::setw(6) << x_ode8 << std::endl;

        x_ode1i = ode1be(x_ode1i, t, sim_options.h, test_system, options_ode1i);
        std::cout << "ode1be    x: " << std::setw(6) << x_ode1i << std::endl;
        std::cout << "number of iterations: ";
        for (const auto& i : options_ode1i.iter_result)
            std::cout << i << " ";
        std::cout << std::endl;
        
        x_ode14x = ode14x(x_ode14x, t, sim_options.h, test_system, options_ode14x);
        std::cout << "ode14x   x: " << std::setw(6) << x_ode14x << std::endl;
        std::cout << "number of iterations: ";
        for (const auto& i : options_ode14x.iter_result)
            std::cout << i << " ";
        std::cout << std::endl;

        std::cout << std::endl;
    }
}
