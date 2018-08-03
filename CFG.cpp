#include <levmar/levmar.h>
#include <cmath>
#include <complex>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>
#include <Faddeeva.hh>

typedef std::vector<double> VD;
//static VD x64 = {0.0243502926634244325089558, 0.0729931217877990394495429,
//                 0.1214628192961205544703765, 0.1696444204239928180373136,
//                 0.2174236437400070841496487, 0.2646871622087674163739642,
//                 0.3113228719902109561575127, 0.3572201583376681159504426,
//                 0.4022701579639916036957668, 0.4463660172534640879849477,
//                 0.4894031457070529574785263, 0.5312794640198945456580139,
//                 0.5718956462026340342838781, 0.6111553551723932502488530,
//                 0.6489654712546573398577612, 0.6852363130542332425635584,
//                 0.7198818501716108268489402, 0.7528199072605318966118638,
//                 0.7839723589433414076102205, 0.8132653151227975597419233,
//                 0.8406292962525803627516915, 0.8659993981540928197607834,
//                 0.8893154459951141058534040, 0.9105221370785028057563807,
//                 0.9295691721319395758214902, 0.9464113748584028160624815,
//                 0.9610087996520537189186141, 0.9733268277899109637418535,
//                 0.9833362538846259569312993, 0.9910133714767443207393824,
//                 0.9963401167719552793469245, 0.9993050417357721394569056};
//static VD w64 = {0.0486909570091397203833654, 0.0485754674415034269347991,
//                 0.0483447622348029571697695, 0.0479993885964583077281262,
//                 0.0475401657148303086622822, 0.0469681828162100173253263,
//                 0.0462847965813144172959532, 0.0454916279274181444797710,
//                 0.0445905581637565630601347, 0.0435837245293234533768279,
//                 0.0424735151236535890073398, 0.0412625632426235286101563,
//                 0.0399537411327203413866569, 0.0385501531786156291289625,
//                 0.0370551285402400460404151, 0.0354722132568823838106931,
//                 0.0338051618371416093915655, 0.0320579283548515535854675,
//                 0.0302346570724024788679741, 0.0283396726142594832275113,
//                 0.0263774697150546586716918, 0.0243527025687108733381776,
//                 0.0222701738083832541592983, 0.0201348231535302093723403,
//                 0.0179517157756973430850453, 0.0157260304760247193219660,
//                 0.0134630478967186425980608, 0.0111681394601311288185905,
//                 0.0088467598263639477230309, 0.0065044579689783628561174,
//                 0.0041470332605624676352875, 0.0017832807216964329472961};





static VD x128 = {0.0122236989606157641980521,0.0366637909687334933302153,0.0610819696041395681037870,0.0854636405045154986364980,0.1097942311276437466729747,0.1340591994611877851175753,0.1582440427142249339974755,0.1823343059853371824103826,0.2063155909020792171540580,0.2301735642266599864109866,0.2538939664226943208556180,0.2774626201779044028062316,0.3008654388776772026671541,0.3240884350244133751832523,0.3471177285976355084261628,0.3699395553498590266165917,0.3925402750332674427356482,0.4149063795522750154922739,0.4370245010371041629370429,0.4588814198335521954490891,0.4804640724041720258582757,0.5017595591361444642896063,0.5227551520511754784539479,0.5434383024128103634441936,0.5637966482266180839144308,0.5838180216287630895500389,0.6034904561585486242035732,0.6228021939105849107615396,0.6417416925623075571535249,0.6602976322726460521059468,0.6784589224477192593677557,0.6962147083695143323850866,0.7135543776835874133438599,0.7304675667419088064717369,0.7469441667970619811698824,0.7629743300440947227797691,0.7785484755064119668504941,0.7936572947621932902433329,0.8082917575079136601196422,0.8224431169556438424645942,0.8361029150609068471168753,0.8492629875779689691636001,0.8619154689395484605906323,0.8740527969580317986954180,0.8856677173453972174082924,0.8967532880491581843864474,0.9073028834017568139214859,0.9173101980809605370364836,0.9267692508789478433346245,0.9356743882779163757831268,0.9440202878302201821211114,0.9518019613412643862177963,0.9590147578536999280989185,0.9656543664319652686458290,0.9717168187471365809043384,0.9771984914639073871653744,0.9820961084357185360247656,0.9864067427245862088712355,0.9901278184917343833379303,0.9932571129002129353034372,0.9957927585349811868641612,0.9977332486255140198821574,0.9990774599773758950119878,0.9998248879471319144736081};
static VD w128 = {0.0244461801962625182113259,0.0244315690978500450548486,0.0244023556338495820932980,0.0243585572646906258532685,0.0243002001679718653234426,0.0242273192228152481200933,0.0241399579890192849977167,0.0240381686810240526375873,0.0239220121367034556724504,0.0237915577810034006387807,0.0236468835844476151436514,0.0234880760165359131530253,0.0233152299940627601224157,0.0231284488243870278792979,0.0229278441436868469204110,0.0227135358502364613097126,0.0224856520327449668718246,0.0222443288937997651046291,0.0219897106684604914341221,0.0217219495380520753752610,0.0214412055392084601371119,0.0211476464682213485370195,0.0208414477807511491135839,0.0205227924869600694322850,0.0201918710421300411806732,0.0198488812328308622199444,0.0194940280587066028230219,0.0191275236099509454865185,0.0187495869405447086509195,0.0183604439373313432212893,0.0179603271850086859401969,0.0175494758271177046487069,0.0171281354231113768306810,0.0166965578015892045890915,0.0162550009097851870516575,0.0158037286593993468589656,0.0153430107688651440859909,0.0148731226021473142523855,0.0143943450041668461768239,0.0139069641329519852442880,0.0134112712886163323144890,0.0129075627392673472204428,0.0123961395439509229688217,0.0118773073727402795758911,0.0113513763240804166932817,0.0108186607395030762476596,0.0102794790158321571332153,0.0097341534150068058635483,0.0091830098716608743344787,0.0086263777986167497049788,0.0080645898904860579729286,0.0074979819256347286876720,0.0069268925668988135634267,0.0063516631617071887872143,0.0057726375428656985893346,0.0051901618326763302050708,0.0046045842567029551182905,0.0040162549837386423131943,0.0034255260409102157743378,0.0028327514714579910952857,0.0022382884309626187436221,0.0016425030186690295387909,0.0010458126793403487793129,0.0004493809602920903763943};

typedef std::complex<double> CD;
CD one(1.0, 0.0), zero(0.0, 0.0), two(2.0, 0.0), i(0.0, 1.0);
const double pi = 4.0 * atan(1.0), lb = 0.0, ub = 200, mid = 0.5 * (ub + lb),
             halfRange = 0.5 * (ub - lb);

struct GLsetting {
    int nGrid;
    VD* abs;
    VD* weights;
};

//static GLsetting gl = {64, &x64, &w64};
static GLsetting gl = {128, &x128, &w128};

struct CFPriceData {
    CD CFPrice;
    CD iu;
    CD interU;
    CD xi;
    CD d;
    CD interD;
    CD A1;
    CD A2;
    CD A;
    CD D;
};

struct CFvolData{
    CD CFvol;
    CD G;
    CD F;
};

struct intVIXdata{
    double VIXint;
    CD rawVIXint;
    double atauBar;
    CD G;
    CD F;
    CD inputU;
    CD iu;
};

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

CFPriceData CFprice(CD u, modelPara p, double tau, double S, double r);
double SPXintegrand(double u, modelPara p, double tau, double K, double S,
                    double r);
VD SPXprice(modelPara p, VD tau, double S, VD K, double r, int n);
CFvolData CFvol(CD u, modelPara p, double tau);
intVIXdata VIXintegrand(CD u, modelPara p, double tau, double K, double tbar);
VD VIXprice(modelPara p, VD &tau, double tbar, VD &K, double r, int n);
VD gradSPXintgrand(double u, modelPara p, double tau, double K, double S, double r);
VD gradientSPXprice(modelPara p, double S, double r, int n, VD &tau, VD &K);
VD gradVIXintegrand(CD u, modelPara p, double tau, double K, double tbar);

//For comparing
void showSPXcallPrices(modelPara mp, VD tarr, double S, VD karr, double r, int n);
void printCFvol(modelPara p, double tau, int n);
void printIntegrandVIXoption(modelPara p, double tau, double K, double tbar, int n);
void printVIXcalls(modelPara p, VD tau, double tbar, VD K, double r, int n);
void printSPXgradient(modelPara p, double S, double r, int n, VD tau, VD K);

int main() {
    modelPara mp = {3.0, 0.1, 0.08, -0.8, 0.25};
    double S = 1.0;
    double r = 0.02;
    double tbar = 30/365.0;
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

    mktPara marP = {S, r, tarr, karr};

    //showSPXcallPrices(mp, tarr, S, karr, r, (int)karr.size());
    printCFvol(mp, tarr[1], gl.nGrid>>1);
    printIntegrandVIXoption(mp, tarr[5], karr[5], tbar, gl.nGrid>>1);
    printVIXcalls(mp, tarr, tbar, karr, r, (int)karr.size());
    //printSPXgradient(mp, S, r, (int)karr.size(), tarr, karr);
}


//Functions for comparing results
void printSPXgradient(modelPara p, double S, double r, int n, VD tau, VD K){
    VD gradients = gradientSPXprice(p, S, r, n, tau, K);
    for (int j = 0; j < 200; j++){
        std::cout << std::setprecision(16);
        std::cout << "gradient" << j << "   " <<  gradients[j] << std::endl;
    }
}


void showSPXcallPrices(modelPara mp, VD tarr, double S, VD karr, double r, int n){

    VD SPXprices = SPXprice(mp, tarr, S, karr, r, n);

    for (int j = 0; j < n; j++) {
        std::cout << std::fixed << std::setprecision(16);
        std::cout << "strike_price  " << karr[j] << std::endl;
        std::cout << "maturity      " << tarr[j] << std::endl;
        std::cout << "SPXcall_price " << SPXprices[j] << std::endl
                  << std::endl;
    }
}

void printVIXcalls(modelPara p, VD tau, double tbar, VD K, double r, int n){
    VD VIXcalls = VIXprice(p, tau, tbar, K, r, n);
    for (int j = 0; j < n; j++) {
        std::cout << std::setprecision(16);
        std::cout << "strike_price  " << K[j] << std::endl;
        std::cout << "maturity      " << tau[j] << std::endl;
        std::cout << "VIXcall_price " << VIXcalls[j] << std::endl
                  << std::endl;
    }
}

void printCFvol(modelPara p, double tau, int n){
    CFvolData CFv;
    VD u = *gl.abs;
    for(int j = 0; j < n; j++){
        std::cout << std::setprecision(16);
        std::cout << "u " << u[j] <<  std::endl;
        CFv = CFvol(u[j]*one, p, tau);
        std::cout << "CFvol " << CFv.CFvol << std::endl;
    }
}

void printIntegrandVIXoption(modelPara p, double tau, double K, double tbar, int n){
    double VIXint;
    VD u = *gl.abs;
    for (int j = 0; j < n; j++){
        std::cout << std::setprecision(16);
        std::cout << "u " << u[j] << ", ";
        VIXint = VIXintegrand(u[j], p, tau, K, tbar).VIXint;
        std::cout << "VIXintegrand " << VIXint << std::endl;
    }
}




//Functions to keep
CFPriceData CFprice(CD u, modelPara p, double tau, double S, double r) {
    double var = pow(p.sigma, 2);

    CD iu = i * u;
    CD interU = pow(u, 2) + iu;
    CD xi = p.k - p.sigma * p.rho * iu;
    CD d = sqrt(pow(xi, 2) + var * interU);
    CD interD = d * tau * 0.5;
    CD A1 = interU * sinh(interD);
    CD A2 = (d * cosh(interD) + xi * sinh(interD)) / p.v0;
    CD A = A1 / A2;
    CD D = log(d / (p.v0 * A2)) + p.k * tau * 0.5;
    double tmp = p.k * p.vbar / p.sigma;

    CD CF = exp(iu * (log(S) + r * tau) - tmp * p.rho * tau * iu - A +
                2 * tmp * D / p.sigma);


    CFPriceData ret = {CF, iu, interU, xi, d, interD, A1, A2, A, D};
    return ret;
}

VD gradSPXintgrand(double u, modelPara p, double tau, double K, double S, double r){
    VD gradSPXint;
    gradSPXint.reserve(5);
    CD inputU = u - i*0.5;
    double var = pow(p.sigma, 2);
    CFPriceData CFP = CFprice(inputU, p, tau, S, r);
    CD B = exp(CFP.D);
    CD tmpBase = p.k*tau*CFP.iu/p.sigma;
    double tmp1 = 2*p.k*p.vbar/var;
    CD dOverA2 = CFP.d/CFP.A2;
    CD interXi = 2.0 + tau*CFP.xi;
    CD revA2 = 1.0/CFP.A2;
    CD AoverA2 = CFP.A / CFP.A2;
    CD coshInterD = cosh(CFP.interD);

    //Partials for rho
    CD d_rho = -CFP.xi*p.sigma*CFP.iu/CFP.d;
    CD A1_rho = -CFP.iu*CFP.interU*tau*CFP.xi*p.sigma/(2.0*CFP.d)*coshInterD;
    CD A2_rho = -p.sigma*CFP.iu*interXi/(2.0*CFP.d*p.v0)
        *(CFP.xi*coshInterD + CFP.d*sinh(CFP.interD));
    CD A_rho = revA2*A1_rho - AoverA2*A2_rho;
    CD B_rho = exp(p.k*tau*0.5)/p.v0*(revA2*d_rho-dOverA2/CFP.A2*A2_rho);

    //Partials for k
    CD B_k = i/(p.sigma*inputU)*B_rho + B*tau*0.5;

    //Partials for sigma
    CD d_sigma = (p.rho/p.sigma - 1.0/CFP.xi)*d_rho+p.sigma*pow(inputU, 2)/CFP.d;
    CD A1_sigma = CFP.interU*tau*d_sigma*coshInterD/2.0;
    CD A2_sigma = p.rho/p.sigma*A2_rho - interXi/(p.v0*tau*CFP.xi*CFP.iu)
        *A1_rho + p.sigma*tau*CFP.A1/(2.0*p.v0);
    CD A_sigma = revA2*A1_sigma - AoverA2*A2_sigma;

    //Gradient of CF price
    CD v0Par = (-CFP.A / p.v0) * CFP.CFPrice;
    CD vbarPar = (tmp1/p.vbar*CFP.D - tmpBase*p.rho) * CFP.CFPrice;
    CD rhoPar = (-A_rho+tmp1/CFP.d*(d_rho-dOverA2*A2_rho)-tmpBase*p.vbar)*CFP.CFPrice;
    CD kPar = (1.0/(p.sigma*CFP.iu)*A_rho + tmp1/p.k*CFP.D + tmp1/B*B_k
        - tmpBase/p.k*p.vbar*p.rho) * CFP.CFPrice;
    CD sigmaPar = (-A_sigma - 2.0*tmp1/p.sigma*CFP.D + tmp1/CFP.d*(d_sigma
        - dOverA2*A2_sigma) + tmpBase*p.vbar*p.rho/p.sigma) * CFP.CFPrice;

    double x = log(S);
    double rT = r * tau;
    double kappa = x - log(K) + rT;
    CD integrand1 = exp(i * u * kappa - i * inputU * (x + rT));
    double integrand2 = pow(u, 2) + 0.25;
    gradSPXint[0] = real(integrand1 * v0Par) / integrand2;
    gradSPXint[1] = real(integrand1 * vbarPar) / integrand2;
    gradSPXint[2] = real(integrand1 * rhoPar) / integrand2;
    gradSPXint[3] = real(integrand1 * kPar) / integrand2;
    gradSPXint[4] = real(integrand1 * sigmaPar) / integrand2;
    return gradSPXint;
}

CFvolData CFvol(CD u, modelPara p, double tau) {
    double var = pow(p.sigma, 2);

    CD iu = u * i;
    double ktauH = p.k*tau*0.5;
    CD G = cosh(ktauH) + (1.0 - var*iu/p.k)*sinh(ktauH);
    CD F = p.v0*iu/G * exp(-ktauH);

    CD CF = pow(exp(ktauH)/G, 2*p.k*p.vbar/var) * exp(F);
    CFvolData ret = {CF, G, F};

    return ret;
}

double SPXintegrand(double u, modelPara p, double tau, double K, double S,
                    double r) {
    CD iu = i * u;
    CD inputU = u - i * 0.5;
    double x = log(S);
    double rT = r * tau;
    double kappa = x - log(K) + rT;
    CD integrand1 = exp(iu * kappa - i * inputU * (x + rT));
    CFPriceData tmp = CFprice(inputU, p, tau, S, r);
    CD CFpri = tmp.CFPrice;
    double SPXint = real(integrand1 * CFpri) / (pow(u, 2) + 0.25);

    return SPXint;
}

intVIXdata VIXintegrand(CD u, modelPara p, double tau, double K, double tbar){
    CD iu = i * u;
    double k = p.k;
    double vbar = p.vbar;
    double atauBar = (1.0 - exp(-tbar*k))/k;
    CD inputU = -u * atauBar/tbar;
    double btauBar = vbar*(tbar - atauBar);
    CFvolData tmp = CFvol(inputU, p, tau);
    CD CFvola = tmp.CFvol;
    CD part1 = exp(-iu*btauBar/tbar);
    CD part2 = 1.0 - Faddeeva::erf(K/100.0 * sqrt(-iu));
    CD part3 = pow(-iu, 3/2.0);
    CD rawInt = CFvola*part1*part2/part3;

    double VIXint = real(rawInt);

    intVIXdata ret = {VIXint, rawInt, atauBar, tmp.G, tmp.F, inputU, iu};

    return ret;
}

VD gradVIXintegrand(CD u, modelPara p, double tau, double K, double tbar){
    intVIXdata inter = VIXintegrand(u, p, tau, K, tbar);
    double var = pow(p.sigma, 2);
    double tmp1 = 2.0*p.k/var;
    double tmp2 = p.k*tau*0.5;
    double tmp3 = exp(tmp2);
    CD tmp4 = pow(inter.G, 2);
    CD iU = i * inter.inputU;
    CD iuOverTbar = inter.iu/tbar;

    //equation (3.34)
    CD G_sigma = -2.0*p.sigma*iU/p.k*sinh(tmp2);

    CD h_v0 = inter.F/p.v0;
    CD h_vbar = tmp1*log(tmp3/inter.G);
    CD h_sigma = -2.0*p.vbar/p.sigma*h_vbar - tmp1*p.vbar/inter.G*G_sigma
        - p.v0*iU/(tmp4*tmp3)*G_sigma;
    CD h_k = -p.sigma/(2.0*p.k)*h_sigma + p.vbar*tau*iU/(inter.G*tmp3)
        - p.v0*inter.inputU*tau/(2.0*p.k*tmp4)*(2.0*p.k*i + inter.inputU*var);

    //equation (3.16)
    double atauBar_k = (tbar - inter.atauBar*(p.k*tbar + 1))/p.k;
    double btauBar_vbar = p.vbar - inter.atauBar;
    double btauBar_k = -p.vbar*atauBar_k;

    //equation (3.23), (3.24)
    CD G_U = (inter.G - tmp3)/inter.inputU;
    CD hPrime = h_k - u/tbar*(inter.F/inter.inputU - 1.0/inter.G*(tmp1*p.vbar
        + inter.F)*G_U)*atauBar_k;

    //equation (3.29)
    CD H_vbar = h_vbar - iuOverTbar*btauBar_vbar;
    CD H_k = hPrime - iuOverTbar*btauBar_k;


    VD ret;
    ret.reserve(5);

    //equation (3.28)
    ret[0] = real(h_v0*inter.rawVIXint);
    ret[1] = real(H_vbar*inter.rawVIXint);
    ret[2] = 0.0;
    ret[3] = real(H_k*inter.rawVIXint);
    ret[4] = real(h_sigma*inter.rawVIXint);


    return ret;
}


VD SPXprice(modelPara p, VD tau, double S, VD K, double r, int n) { //tau and K can be passed by reference, but we will see, may be I'll make them const or static
    int nGrid = gl.nGrid;

    double up_u, down_u, upInt, downInt, strike, T, rT, glCollect,
        glInt, SPXcall;
    nGrid = nGrid >> 1;
    VD u = *gl.abs;
    VD w = *gl.weights;

    VD SPXs;
    SPXs.reserve(n);

    for (int j = 0; j < n; j++) {
        strike = K[j];
        T = tau[j];
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

    return SPXs; //Here can return an adress of the VD.
}

VD gradientSPXprice(modelPara p, double S, double r, int n, VD &tau, VD& K){
    int nGrid = gl.nGrid>>1;
    VD u = *gl.abs;
    VD w = *gl.weights;

    VD gradSPX, upInt, downInt, glCollect;
    gradSPX.reserve(5*n);
    upInt.reserve(5);
    downInt.reserve(5);
    glCollect.reserve(5);

    int k, l;
    double strike, T, rT, up_u, down_u;
    for(k=l=0; l < n; l++){
        strike = K[l];//do I need to use *K?
        T = tau[l];
        rT = r*T;
        for (int cc = 0; cc < 5; cc++){
            glCollect[cc] = 0.0;
        }
        for (int count=0; count < nGrid; count++){
            up_u = mid + u[count] * halfRange;
            down_u = mid - u[count] * halfRange;
            upInt = gradSPXintgrand(up_u, p, T, strike, S, r);
            downInt = gradSPXintgrand(down_u, p, T, strike, S, r);
            for (int j = 0; j < 5; j++)
                glCollect[j] += w[count]*(upInt[j] + downInt[j]);
        }
        for (int p = 0; p < 5; p++){
            glCollect[p] = glCollect[p]*halfRange;
            gradSPX[k++] = - sqrt(strike * S) * exp(-rT * 0.5) / pi * glCollect[p];
        }
    }

    return gradSPX;
}

//For this pricing, the u is not a real but can be a complex. See smile page
//7, equation (11). Currently, u is real which means im(u) is 0, but in the
//paper it said it should be > 0. Do I choose one im(u) and how to choose?
//So does this numerical integration still work? And what should the upper
//boundary be?
VD VIXprice(modelPara p, VD &tau, double tbar, VD &K, double r, int n){
    int nGrid = gl.nGrid;

    double up_u, down_u, upInt, downInt, strike, T, discount, glCollect,
        glInt, VIXcall;
    nGrid = nGrid >> 1;
    VD u = *gl.abs;
    VD w = *gl.weights;

    VD VIXs;
    VIXs.reserve(n);

    for (int j = 0; j < n; j++){
        strike = K[j];
        T = tau[j];
        discount = exp(-r*T);
        glCollect = 0.0;
        for (int count = 0; count < nGrid; count++){
            up_u = mid + u[count] * halfRange;
            down_u = mid - u[count] * halfRange;
            upInt = VIXintegrand(up_u+i, p, T, strike, tbar).VIXint; //What should be the complex part be? By adding this magic number I already made the price positive, but how to choose the complex part?
            downInt = VIXintegrand(down_u+i, p, T, strike, tbar).VIXint;
            glCollect += w[count] * (upInt + downInt);
        }

        glInt = halfRange * glCollect;
        VIXcall = 50.0 * discount/sqrt(pi) * glInt;
        VIXs[j] = VIXcall;
    }
    return VIXs;
}









