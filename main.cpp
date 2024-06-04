#include <iostream>
#include <cmath>
#include <limits>
#include <stdio.h>

class DoubleDouble {
private:
    double hi; // Most significant part
    double lo; // Least significant part

public:
    DoubleDouble(double h = 0.0, double l = 0.0) : hi(h), lo(l) {}

    // Quick two-sum operation
    static DoubleDouble quick_two_sum(double a, double b) {
        double s = a + b;
        double e = b - (s - a);
        return DoubleDouble(s, e);
    }

    // Two-sum operation
    static DoubleDouble two_sum(double a, double b) {
        double s = a + b;
        double v = s - a;
        double e = (a - (s - v)) + (b - v);
        return DoubleDouble(s, e);
    }

    // Two-product operation
    static DoubleDouble two_prod(double a, double b) {
        double p = a * b;
        double e = std::fma(a, b, -p);
        return DoubleDouble(p, e);
    }

    // Addition operation
    static DoubleDouble add(DoubleDouble a, DoubleDouble b) {
        DoubleDouble s = two_sum(a.hi, b.hi);
        s.lo += (a.lo + b.lo);
        return quick_two_sum(s.hi, s.lo);
    }

    // Multiplication operation
    static DoubleDouble multiply(DoubleDouble a, DoubleDouble b) {
        DoubleDouble p1 = two_prod(a.hi, b.hi);
        p1.lo += a.hi * b.lo + a.lo * b.hi;
        return quick_two_sum(p1.hi, p1.lo);
    }

    // Division operation
    static DoubleDouble divide(DoubleDouble a, DoubleDouble b) {
        double q1 = a.hi / b.hi;
        DoubleDouble p = multiply(DoubleDouble(q1), b);
        DoubleDouble r = add(a, DoubleDouble(-p.hi, -p.lo));
        double q2 = r.hi / b.hi;
        return add(DoubleDouble(q1), DoubleDouble(q2));
    }

    // Exponential function
    static DoubleDouble exp(DoubleDouble x) {

        if (x.hi == 0.0 && x.lo == 0.0) {
            return DoubleDouble(1.0, 0.0);
        }

        bool negative = x.hi < 0.0;
        if (negative) {
            x.hi = -x.hi;
            x.lo = -x.lo;
        }

        DoubleDouble sum(1.0, 0.0);
        DoubleDouble term(1.0, 0.0);
        int n = 1;

        while (true) {
            term = divide(multiply(term, x), DoubleDouble(n, 0.0));
            DoubleDouble newSum = add(sum, term);

            if (std::isinf(newSum.hi) || std::isnan(newSum.hi)) {
                if (negative) {
                    return DoubleDouble(0.0, 0.0);
                } else {
                    return DoubleDouble(std::numeric_limits<double>::infinity(), 0.0);
                }
            }

            if (newSum.hi == sum.hi && newSum.lo == sum.lo) {
                break;
            }
            sum = newSum;
            ++n;
        }

        if (negative) {
            return divide(DoubleDouble(1.0, 0.0), sum);
        }

        return sum;
    }

    // Getter for the hi part
    double get_hi() const {
        return hi;
    }

    // Getter for the lo part
    double get_lo() const {
        return lo;
    }
};

int main() {
    DoubleDouble x_values[] = {DoubleDouble(0.0, 0.0), DoubleDouble(1.0, 0.0), DoubleDouble(2.0, 0.0),
                               DoubleDouble(3.0, 0.0), DoubleDouble(4.0, 0.0), DoubleDouble(5.0, 0.0),
                               DoubleDouble(6.0, 0.0), DoubleDouble(-650.0, 0.0), DoubleDouble(700.0, 0.0),
                               DoubleDouble(-1000, 0.0), DoubleDouble(1000, 0.0)};

    for (int i = 0; i < 11; ++i) {
        DoubleDouble result = DoubleDouble::exp(x_values[i]);
        std::cout << "exp(" << x_values[i].get_hi() << ") = ";
        printf("hi: %.40e   lo: %.40e\n", result.get_hi(), result.get_lo());
    }

    return 0;
}
//exp(-650) from wolfram 3.9754497359086468077890997537948254523324502696238e-31 3.975449735908647e-31 -3.3648383059169216e-48 + 
//exp(700) from wolfram 1.0142320547350045094553295952312676152046795722431e+304 1.0142320547350045e+304 1.6666571920734673e+287 +