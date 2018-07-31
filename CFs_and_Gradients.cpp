#include <levmar/levmar.h>
#include <cmath>
#include <complex>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

typedef std::vector<double> VD;
static VD x64 = {0.0243502926634244325089558, 0.0729931217877990394495429,
                 0.1214628192961205544703765, 0.1696444204239928180373136,
                 0.2174236437400070841496487, 0.2646871622087674163739642,
                 0.3113228719902109561575127, 0.3572201583376681159504426,
                 0.4022701579639916036957668, 0.4463660172534640879849477,
                 0.4894031457070529574785263, 0.5312794640198945456580139,
                 0.5718956462026340342838781, 0.6111553551723932502488530,
                 0.6489654712546573398577612, 0.6852363130542332425635584,
                 0.7198818501716108268489402, 0.7528199072605318966118638,
                 0.7839723589433414076102205, 0.8132653151227975597419233,
                 0.8406292962525803627516915, 0.8659993981540928197607834,
                 0.8893154459951141058534040, 0.9105221370785028057563807,
                 0.9295691721319395758214902, 0.9464113748584028160624815,
                 0.9610087996520537189186141, 0.9733268277899109637418535,
                 0.9833362538846259569312993, 0.9910133714767443207393824,
                 0.9963401167719552793469245, 0.9993050417357721394569056};
static VD w64 = {0.0486909570091397203833654, 0.0485754674415034269347991,
                 0.0483447622348029571697695, 0.0479993885964583077281262,
                 0.0475401657148303086622822, 0.0469681828162100173253263,
                 0.0462847965813144172959532, 0.0454916279274181444797710,
                 0.0445905581637565630601347, 0.0435837245293234533768279,
                 0.0424735151236535890073398, 0.0412625632426235286101563,
                 0.0399537411327203413866569, 0.0385501531786156291289625,
                 0.0370551285402400460404151, 0.0354722132568823838106931,
                 0.0338051618371416093915655, 0.0320579283548515535854675,
                 0.0302346570724024788679741, 0.0283396726142594832275113,
                 0.0263774697150546586716918, 0.0243527025687108733381776,
                 0.0222701738083832541592983, 0.0201348231535302093723403,
                 0.0179517157756973430850453, 0.0157260304760247193219660,
                 0.0134630478967186425980608, 0.0111681394601311288185905,
                 0.0088467598263639477230309, 0.0065044579689783628561174,
                 0.0041470332605624676352875, 0.0017832807216964329472961};

typedef std::complex<double> CD;
CD one(1.0, 0.0), zero(0.0, 0.0), two(2.0, 0.0), i(0.0, 1.0);
const double pi = 4.0 * atan(1.0), lb = 0.0, ub = 200, mid = 0.5 * (ub + lb),
             halfRange = 0.5 * (ub - lb);

struct GLsetting {
    int nGrid;
    VD* abs;
    VD* weights;
};

static GLsetting gl = {64, &x64, &w64};

struct mktPara {
    double S;
    double r;
    VD T;
    VD K;
};

struct modelPara {
    double k;
    double vbar;
    double v0;
    double rho;
    double sigma;
};

CD CFprice(CD u, modelPara p, double tau, double S, double r);
double SPXintegrand(double u, modelPara p, double tau, double K, double S,
                    double r);
VD SPXprice(modelPara p, VD tau, double S, VD K, double r, int n);

int main() {
    modelPara mp = {3.0, 0.1, 0.08, -0.8, 0.25};
    double S = 1.0;
    double r = 0.02;
    const VD karr = {0.9371, 0.8603, 0.8112, 0.7760, 0.7470, 0.7216, 0.6699,
                     0.6137, 0.9956, 0.9868, 0.9728, 0.9588, 0.9464, 0.9358,
                     0.9175, 0.9025, 1.0427, 1.0463, 1.0499, 1.0530, 1.0562,
                     1.0593, 1.0663, 1.0766, 1.2287, 1.2399, 1.2485, 1.2659,
                     1.2646, 1.2715, 1.2859, 1.3046, 1.3939, 1.4102, 1.4291,
                     1.4456, 1.4603, 1.4736, 1.5005, 1.5328};

    const VD tarr = {0.119047619047619, 0.238095238095238, 0.357142857142857,
                     0.476190476190476, 0.595238095238095, 0.714285714285714,
                     1.07142857142857,  1.42857142857143,  0.119047619047619,
                     0.238095238095238, 0.357142857142857, 0.476190476190476,
                     0.595238095238095, 0.714285714285714, 1.07142857142857,
                     1.42857142857143,  0.119047619047619, 0.238095238095238,
                     0.357142857142857, 0.476190476190476, 0.595238095238095,
                     0.714285714285714, 1.07142857142857,  1.42857142857143,
                     0.119047619047619, 0.238095238095238, 0.357142857142857,
                     0.476190476190476, 0.595238095238095, 0.714285714285714,
                     1.07142857142857,  1.42857142857143,  0.119047619047619,
                     0.238095238095238, 0.357142857142857, 0.476190476190476,
                     0.595238095238095, 0.714285714285714, 1.07142857142857,
                     1.42857142857143};

    mktPara mP = {S, r, tarr, karr};

    VD SPXprices = SPXprice(mp, tarr, S, karr, r, karr.size());

    for (int j = 0; j < karr.size(); j++) {
        std::cout << "strike price =  " << karr[j] << std::endl;
        std::cout << "maturity =      " << tarr[j] << std::endl;
        std::cout << "SPXcall price = " << SPXprices[j] << std::endl
                  << std::endl;
    }
}

CD CFprice(CD u, modelPara p, double tau, double S, double r) {
    double k = p.k;
    double vbar = p.vbar;
    double v0 = p.v0;
    double rho = p.rho;
    double sigma = p.sigma;

    double var = pow(sigma, 2);

    CD iu = i * u;
    CD interU = pow(u, 2) + iu;
    CD xi = k - sigma * rho * iu;
    CD d = sqrt(pow(xi, 2) + var * interU);
    CD interD = d * tau * 0.5;
    CD A1 = interU * sinh(interD);
    CD A2 = (d * cosh(interD) + xi * sinh(interD)) / v0;
    CD A = A1 / A2;
    CD D = log(d / (v0 * A2)) + k * tau * 0.5;
    double tmp = k * vbar / sigma;

    CD CF = exp(iu * (log(S) + r * tau) - tmp * rho * tau * iu - A +
                2 * tmp * D / sigma);

    return CF;
}

double SPXintegrand(double u, modelPara p, double tau, double K, double S,
                    double r) {
    CD iu = i * u;
    CD interU = u - i * 0.5;
    double x = log(S);
    double rT = r * tau;
    double kappa = log(S / K) + rT;
    CD integrand1 = exp(iu * kappa - i * interU * (x + rT));
    CD CFpri = CFprice(u * one, p, tau, S, r);
    double SPXint = real(integrand1 * CFpri) / (pow(u, 2) + 0.25);

    return SPXint;
}

VD SPXprice(modelPara p, VD tau, double S, VD K, double r, int n) {
    int nGrid = gl.nGrid;

    double up_u, down_u, upInt, downInt, strike, T, discountF, rT, glCollect,
        glInt, SPXcall;
    nGrid = nGrid >> 1;
    VD u = *gl.abs;
    VD w = *gl.weights;

    VD SPXs;
    SPXs.reserve(n);

    for (int j = 0; j < n; j++) {
        strike = K[j];
        T = tau[j];
        discountF = exp(-r * T);
        rT = r * T;
        glCollect = 0.0;

        for (int count = 0; count < nGrid; count++) {
            up_u = mid + u[count] * halfRange;
            down_u = mid - u[count] * halfRange;
            upInt = SPXintegrand(up_u, p, T, strike, S, r);
            downInt = SPXintegrand(down_u, p, T, strike, S, r);
            glCollect += w[count] * (upInt + downInt);
        }

        glInt = halfRange * glCollect;
        SPXcall = S - sqrt(strike * S) * exp(-rT * 0.5) / pi * glInt;
        SPXs[j] = SPXcall;
    }

    return SPXs;
}
