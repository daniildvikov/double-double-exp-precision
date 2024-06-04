#include <iostream>
#include <cmath>
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
    static DoubleDouble exp_dd(const DoubleDouble& x) {
        const double ln2_hi = 6.93147180369123816490e-01;
        const double ln2_lo = 1.90821492927058770002e-10;
        const double inv_ln2 = 1.44269504088896338700e+00;

        // Decompose x into n and r so that x = n * ln(2) + r
        int n = static_cast<int>(x.hi * inv_ln2 + (x.hi >= 0 ? 0.5 : -0.5));
        DoubleDouble r = add(x, DoubleDouble(-n * ln2_hi, -n * ln2_lo));

        // Compute exp(r) using Taylor series
        DoubleDouble result(1.0, 0.0);
        DoubleDouble term(1.0, 0.0);    
        DoubleDouble r2 = r;

        for (int i = 1; i < 100; ++i) {
            term = multiply(term, r2);
            term.hi /= i;
            result = add(result, term);
            if (std::abs(term.hi) < std::abs(result.hi) * std::numeric_limits<double>::epsilon()) {
                break;
            }
        }

        // Scale the result by 2^n
        if (n > 0) {
            for (int i = 0; i < n; ++i) {
                result = multiply(result, DoubleDouble(2.0, 0.0));
            }
        } else {
            for (int i = 0; i < -n; ++i) {
                result = multiply(result, DoubleDouble(0.5, 0.0));
            }
        }

        return result;
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
    DoubleDouble x(1.0, 0.0);
    DoubleDouble result = DoubleDouble::exp_dd(x);
    printf("exp(%.0f) = hi : %.40e lo: %.40e\n", x.get_hi(), result.get_hi(), result.get_lo());
    return 0;
}
//exp(1) = hi: 2.7182818284590450907955982984276488423347e+00   lo: 1.4456468917292496647933303150418560929681e-16