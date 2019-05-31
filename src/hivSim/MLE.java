package hivSim;

import java.math.BigDecimal;
import java.math.BigInteger;
import java.math.MathContext;
import java.math.RoundingMode;
import java.util.*;

import static hivSim.CSV.rValues;
import static java.lang.Math.*;
import static nr.ran.NRUtil.*;

/**
 * Contains methods related to calculating the maximum likelihood estimate and
 * standard deviations of the basic reproduction number. Contains a global
 * class variable 'tolerance' which controls the number of significant figures
 * desired; the default value is 9. The method setTolerance() in CSV
 * should be used to set this value in the simulation.
 */

public class MLE {
    private static int tolerance;
    private static double[][] Binom;
    private static double[][][] Skmn;
    //private static MathContext mc;

    //public static void setMC(){
    //    mc = new MathContext(tolerance, RoundingMode.HALF_DOWN);
    //}

    /**
     * Generates a lookup table with default values equal to -1.
     *
     * @param n the maximum index to be looked up.
     */
    public static void populateBinom(int n){
        Binom = new double[n+1][n+1];
        for (int i=0; i<=n; i++){
            for (int j = 0; j <= n; j++)
                Binom[i][j] = -1;
        }
    }

    /**
     * Sets the class variable 'tolerance' which controls the number of
     * significant digits used in calculations. Note this method should not
     * be called except by CSV.setTolerance() which globally controls the
     * tolerance.
     *
     * @param i the number of significant figures desired.
     */
    public static void setTolerance(int i){
        tolerance=i;
    }

    /**
     * Directly estimates r0 using the Euler-Maclaurin approximation in the
     * maximum likelihood estimation. Does not maximize a likelihood from an
     * array and instead directly calculates a numerical value corresponding to
     * the maximum likelihood.
     *
     * @param M the number of Viruses sampled.
     * @param mu the expected number of mutations per generation.
     * @param etaTilda the folded site frequency spectrum.
     * @return a double representing the estimate of r0.
     * @see <img src="./doc-files/rHat_MLE_EM.png"/>
     */
    public static double rHatEM(double M, double mu, List<Integer> etaTilda){
        double sum = 0;
        for (int m = 2; m <= floor(M/2.0); m++){
            sum = sum + etaTilda.get(m-1); //arrays start from index 0 instead of 1
        }
        double term1 = mu*M*(1 - (1/(M-2)));
        return exp(term1/sum);
    }

    /**
     * Calculates the estimated standard deviation in the natural log of the
     * estimated r0 using only the estimated r0 and the folded site-frequency
     * spectrum.
     *
     * @param rHat the estimated r0.
     * @param M the number of Viruses in the sample.
     * @param etaTilda the folded site-frequency spectrum.
     * @return the estimated standard deviation in the natural log of the
     *         estimated r0.
     * @see <img src="./doc-files/sHat_MLE_EM.png"/>
     */
    public static double sigmaHatEM(double rHat, int M, List<Integer> etaTilda){
        double log = log(rHat);
        double sum = 0;
        for (int m = 2; m <= floor(M/2.0); m++){
            sum = sum + etaTilda.get(m-1); //arrays start from 0
        }
        double sqrt = sqrt(sum);
        return log/sqrt;
    }

    /**
     * Calculates the maximum likelihood estimate for the basic reproduction
     * number.
     *
     * @param etaTilda the folded site-frequency spectrum.
     * @param M the number of Viruses in the sample.
     * @param mu empirical number of mutations per generation.
     * @param aTildas switch controlling whether aTildas or aDeltas are used in
     *        likelihood calculations, with true corresponding to aTilda.
     * @return the maximum likelihood estimate for r0.
     */
    public static double rHat(List<Integer> etaTilda, int M, double mu, boolean aTildas){
        List<Double> rValues = rValues(1, 10, 0.1);
        double guess = gridSearch(rValues, mu, M, etaTilda, true);
        //List<Double> bounds = goldenBounds(1.05, guess*4, mu, M, etaTilda, aTildas);
        return goldenSection(1.05, guess, guess*4, M, mu, etaTilda, aTildas);
    }

    /**
     * Calculates a double value representing an estimate of the standard
     * deviation in ln(rHat) for the maximum likelihood estimate of r0.
     *
     * @param rHat the maximum likelihood estimate of r0.
     * @param M the sample size.
     * @param mu empirical number of mutations per generation.
     * @param etaTilda folded site frequency spectrum.
     * @param aTildas switch controlling whether aTildas or aDeltas are used in
     *        likelihood calculations, with true corresponding to aTilda.
     * @return an estimate of the standard deviation in ln(rHat) for the
     *         maximum likelihood estimation of r0.
     * @see <img src="./doc-files/FisherInformation.png"/>
     * @see <img src="./doc-files/varLnRHat.png"/>
     */

    public static double sigmaHat(double rHat, int M, double mu, List<Integer> etaTilda, boolean aTildas){
        double term, sum = 0;
        if (aTildas){
            for (int m = 2; m <= floor(M/2.0); m++){
                double muA = mu*aTilda2Prime(rHat, M, m);
                double eta = etaTilda.get(m-1); //index starts at 0
                double ddr =     aTildaPrime(rHat, M, m)
                        /aTilda(rHat, M, m);
                double ddr2 = pow(ddr, 2);
                double d2dr2 =  aTilda2Prime(rHat, M, m)
                        /aTilda(rHat, M, m);
                double etaTerm = eta*(ddr2 - d2dr2);
                term = muA + etaTerm;
                sum = sum + term;
            }
        }
        else {
            for (int m = 2; m <= floor(M/2.0); m++) {
                double muA = mu * aDelta2Prime(rHat, M, m);
                double eta = etaTilda.get(m - 1); //index starts at 0
                double ddr = aDeltaPrime(rHat, M, m)
                        / aDelta(rHat, M, m);
                double ddr2 = pow(ddr, 2);
                double d2dr2 = aDelta2Prime(rHat, M, m)
                        / aDelta(rHat, M, m);
                double etaTerm = eta * (ddr2 - d2dr2);
                term = muA + etaTerm;
                sum = sum + term;
            }
        }
        double inverseR = pow(rHat, -1);
        double inverseFisher2 = pow(sum, -2);
        double var = inverseR*inverseFisher2;
        return Math.sqrt(var);
    }

    /**
     * Calculates the log likelihood of a phylogeny producing the given folded
     * site-frequency spectrum, as a function of r.
     *
     * @param mu the number of mutations per generation.
     * @param r the basic reproduction number.
     * @param M the number of Viruses in the sample.
     * @param etaTilda the folded site frequency spectrum.
     * @param aTildas switch controlling whether aTildas or aDeltas are used in
     *        likelihood calculations, with true corresponding to aTilda.
     * @return the value of the log likelihood function at a given basic
     *         reproduction number.
     * @see <img src="./doc-files/lnL(r).png"/>
     */
    public static double logLikelihood(double mu, double r, int M, List<Integer> etaTilda, boolean aTildas){
        double sum = 0;
        if (aTildas){
            for (int m = 2; m <= floor(M/2.0); m++){
                double aTilda_m = aTilda(r, M, m);
                System.out.printf("%n A(%f, %d %d) = %f", r, M, m, aTilda_m);
                sum += -mu*aTilda_m + etaTilda.get(m-1) * log(aTilda_m); //eta uses m-1 because eta[0] is eta_1
            }
        }
        else{
            for (int m = 2; m <= floor(M/2.0); m++){
                double aDelta_m = aDelta(r, M, m);
                sum += -mu*aDelta_m + etaTilda.get(m-1) * log(aDelta_m); //eta uses m-1 because eta[0] is eta_1
            }
        }

        return sum;
    }

    /**
     * Uses a golden-section search to find the maximum likelihood estimate of
     * the basic reproduction number.
     *
     * @param ax initial lower bound.
     * @param bx initial midpoint.
     * @param cx initial upper bound.
     * @param M sample size.
     * @param mu mean number of mutations per generation.
     * @param etaTilda folded site frequency spectrum.
     * @param aTildas switch controlling whether aTildas or aDeltas are used in
     *        likelihood calculations, with true corresponding to aTilda.
     * @return the maximum likelihood estimate for the basic reproduction
     *         number as estimated by golden section minimization.
     * @see <img src="./doc-files/goldenSection.png"/>
     */
    public static double goldenSection(double ax, double bx, double cx, int M, double mu, List<Integer> etaTilda, boolean aTildas)  {
        final double R=0.61803399,C=1.0-R;
        double xmin;
        double tol = pow(10, -tolerance);

        double x1,x2;
        double x0=ax;
        double x3=cx;
        if (abs(cx - bx) > abs(bx - ax)) {
            x1=bx;
            x2=bx+C*(cx - bx);
        } else {
            x2=bx;
            x1=bx-C*(bx-ax);
        }
        double f1= -logLikelihood(mu, x1, M, etaTilda, aTildas);
        double f2= -logLikelihood(mu, x2, M, etaTilda, aTildas);
        while (abs(x3-x0) > tol*(abs(x1)+abs(x2))) {
            if (f2 < f1) {
                double dum = R*x2+C*x3;
                x0=x1; x1=x2; x2=dum;
                f1=f2;
                f2 = -logLikelihood(mu, x2, M, etaTilda, aTildas);
            } else {
                double dum = R*x1+C*x0;
                x3=x2; x2=x1; x1=dum;
                f2=f1;
                f1= -logLikelihood(mu, x1, M, etaTilda, aTildas);
            }
        }
        if (f1 < f2){xmin=x1;}
        else        {xmin=x2;}
        return xmin;
    }

    /**
     * Uses a grid to identify the maximum likelihood estimate of the basic
     * reproduction number.
     *
     * @param rValues the values in the grid to check.
     * @param mu the average number of mutations per generation.
     * @param M the number of Viruses in the sample.
     * @param etaTilda the folded site frequency spectrum.
     * @param aTildas switch controlling whether aTildas or aDeltas are used in
     *        likelihood calculations, with true corresponding to aTilda.
     * @return the value in the input grid which produces the largest (log)
     *         likelihood, i.e. the maximum likelihood estimate for r0.
     */
    public static double gridSearch(List<Double> rValues, double mu, int M, List<Integer> etaTilda, boolean aTildas){
        List<Double> grid = new ArrayList<>();
        for (double r0: rValues){
            grid.add(logLikelihood(mu, r0, M, etaTilda, aTildas));
        }
        int largest = 0;
        for (int i=1; i<grid.size(); i++){
            if (grid.get(i) > grid.get(largest)){largest = i;}
        }
        return rValues.get(largest);
    }

    /**
     * Calculates the (approximate) expected value of a specific m in the
     * ancestor frequency spectrum as a function of the basic reproduction
     * number and sample size.
     *
     * @param r the basic reproduction number.
     * @param M the sample size.
     * @param m the index of the ancestor frequency spectrum to calculate.
     *        Note element [m] is the count of ancestors with m descendants
     *        in the sample. m > 0.
     *
     * @return the expected value of index m in the ancestor frequency spectrum.
     * @see <img src="./doc-files/aDelta1.png"/>
     * @see <img src="./doc-files/aDelta.png"/>
     */
    public static double aDelta(double r, int M, int m) {
        if (m <= 0) {
            throw new IllegalArgumentException("m must be greater than 0.");
        }
        if (m == 1) { //handle m=1
            double sum = 0;
            double term = 1;
            double tol = pow(10, -tolerance);
            int g = 1;
            while (abs(term) > tol * abs(sum)) {
                double rg = pow(r, -g);
                double base = 1 - (1 - rg);
                term = pow(base, M - 1);
                sum = sum + term;
                g++;
            }
            return M * sum;
        }
        return binomial(M, m) * S_kmn(r, 0, m - 1, M - m); //all other m
    }

    /**
     * Calculates the "folded" version of the expected value of the sample
     * ancestor frequency spectrum by taking sums and differences of aDelta.
     *
     * @param r the basic reproduction number.
     * @param M the number of samples.
     * @param m the index of the (folded) ancestor frequency spectrum to
     *        calculate. 0 < m <= floor(M/2).
     * @return the folded expected value of the ancestor frequency spectrum at
     *         a given index m.
     * @see <img src="./doc-files/aTilda.png"/>
     * @see <img src="./doc-files/aTilda1.png"/>
     */
    public static double aTilda(double r, int M, int m){
        if (m > floor(M)/2.0){throw new IllegalArgumentException("0 < m <= floor(M/2)");}
        if (M%2 == 0 & m == floor(M)/2.0){return aDelta(r, M, m);}
        if (m==1){
            return aDelta(r, M, 1) - aDelta(r, M, M-1);
        }
        return aDelta(r, M, m) + aDelta(r, M, M-m);
    }

    /**
     * Calculates the first derivative (with respect to r) of the expected
     * value of the ancestor frequency spectrum at a given reproduction number
     * and index m.
     *
     * @param r the basic reproduction number.
     * @param M the sample size.
     * @param m the index of the ancestor frequency spectrum whose derivative
     *        is required.
     * @return the r-derivative of the expected value of the ancestor frequency
     *         at the given index m and r0.
     * @see <img src="./doc-files/aDeltaPrime.png"/>
     * @see <img src="./doc-files/aDeltaPrime1.png"/>
     */
    public static double aDeltaPrime(double r, int M, int m){
        if (m == 1){
            return -M*(M-1)*(1/r)*S_kmn(r, 1, 1,   M-2);
        }
        double term1 = binomial(M, m)*(1/r);
        double sum1 =  (M-1)*    S_kmn(r, 1, m,   M-m-1);
        double sum2 = -(m-1)*    S_kmn(r, 1, m-1, M-m-1);
        return term1*(sum1+sum2);
    }

    /**
     * Calculates the second derivative (with respect to r) of the expected
     * value of the ancestor frequency spectrum at a given reproduction number
     * and index m.
     *
     * @param r the basic reproduction number.
     * @param M the sample size.
     * @param m the index of the ancestor frequency spectrum whose 2nd
     *        derivative is required.
     * @return the 2nd r-derivative of the expected value of the ancestor
     *         frequency spectrum at the given index m and r0.
     * @see <img src="./doc-files/aDelta2Prime.png"/>
     * @see <img src="./doc-files/aDelta2Prime1.png"/>
     */
    public static double aDelta2Prime(double r, int M, int m){
        if (m == 1){
            double term1 = -M*(M-1)*pow(r, -2);
            double sum1 = (M-1)*S_kmn(r, 2, 2, M-3);
            double sum2 =      -S_kmn(r, 1, 2, M-3);
            double sum3 =      -S_kmn(r, 2, 1, M-3);
            double sum4 =      -S_kmn(r, 1, 1, M-3);
            return term1*(sum1 + sum2 + sum3 + sum4);
        }
        double term1 = binomial(M, m)/(pow(r, 2));
        double sum1 =  (pow(M-1, 2))         *S_kmn(r, 2, m+1, M-m-2);
        double sum2 =  (M-1)                 *S_kmn(r, 1, m+1, M-m-2);
        double sum3 =  (2*m*M - 3*m - M + 2) *S_kmn(r, 2, m,   M-m-2);
        double sum4 = -(m+M-2)               *S_kmn(r, 1, m,   M-m-2);
        double sum5 =  (pow(m-1, 2))         *S_kmn(r, 2, m-1, M-m-2);
        double sum6 =  (m-1)                 *S_kmn(r, 1, m-1, M-m-2);
        return term1*(sum1 + sum2 + sum3 + sum4 + sum5 + sum6);
    }

    /**
     * Calculates the "folded" version of the first derivative of the expected
     * value of the sample ancestor frequency spectrum at a particular basic
     * reproduction number and index m by taking sums and differences of the
     * unfolded derivatives (e.g. aDeltaPrime).
     *
     * @param r the basic reproduction number.
     * @param M the sample size.
     * @param m the index of the ancestor frequency spectrum whose (folded) 1st
     *        derivative is required.
     * @return the folded 1st r derivative of the expected value of the
     *         ancestor frequency spectrum at a given index m and r0.
     * @see <img src="./doc-files/aTildaPrime.png"/>
     * @see <img src="./doc-files/aTildaPrime1.png"/>
     */
    public static double aTildaPrime(double r, int M, int m){
        if (m==1){
            return aDeltaPrime(r, M, m) - aDeltaPrime(r, M, M-m);
        }
        return aDeltaPrime(r, M, m) + aDeltaPrime(r, M, M-m);
    }

    /**
     * Calculates the "folded" version of the second derivative of the
     * expected value of the sample ancestor frequency spectrum at a particular
     * basic reproduction number and index m by taking sums and differences of
     * the unfolded derivatives (e.g. aDelta2Prime).
     *
     * @param r the basic reproduction number.
     * @param M the sample size.
     * @param m the index of the ancestor frequency spectrum whose (folded) 2nd
     *        derivative is required.
     * @return the folded 2nd r-derivative of the expected value of the
     *         ancestor frequency spectrum at a given index m and r0.
     * @see <img src="./doc-files/aTilda2Prime.png"/>
     * @see <img src="./doc-files/aTilda2Prime1.png"/>
     */
    public static double aTilda2Prime(double r, int M, int m){
        if (m==1){
            return aDelta2Prime(r, M, m) - aDelta2Prime(r, M, M-m);
        }
        return aDeltaPrime(r, M, m) + aDelta2Prime(r, M, M-m);
    }

    /**
     * Calculates the sum S_kmn(r) explicitly. The sum is truncated according
     * to the global value tolerance.
     *
     * @param r the basic reproduction number.
     * @param k the power of g in each term of the sum.
     * @param m the power of r^-g in each term of the sum.
     * @param n the power of 1 - r^-g in each term of the sum.
     * @return the sum S_kmn(r) accurate to 9 decimal places.
     * @see <img src="./doc-files/Skmn.png"/>
     */
    public static double S_kmn(double r, int k, int m, int n){
        //if (Skmn[k][m][n] == -1){
            double sum = 0;
            double term = 1;
            int g = 1;
            double tol = pow(10, -tolerance);
            while (abs(term) > tol*abs(sum)) {
                double rg = pow(r, -g);
                term = (pow(g, k))*(pow(rg, m))*(pow(1-rg, n));
                sum = sum + term;
                g++;
            }
            //Skmn[k][m][n] = sum;
            return sum;
        //}
        //return Skmn[k][m][n];
    }

    /**
     * Calculates the coefficient of the x^k term in the polynomial expansion of
     * (1 + x)^n. Stores these values in a [][] to maximize computational
     * efficiency. Note that the output will be in scientific notation with
     * the number of significant digits controlled by tolerance.
     *
     * @param N the binomial power, as in (1 + x)^N.
     * @param K the power of the term whose coefficient is to be calculated, as in x^K.
     * @return the binomial coefficient N choose K.
     */

    public static double binomial(final int N, final int K) {
        BigInteger result = BigInteger.ONE;
        if (Binom[N][K] != -1){
            return Binom[N][K];
        }
        for (int k = 0; k < K; k++) {
            result = result.multiply(BigInteger.valueOf(N-k));
            result = result.divide(  BigInteger.valueOf(k+1));
        }
        Binom[N][K] = result.doubleValue();
        return result.doubleValue();
    }

    //Below this point are functions which are not used in the main body of the
    //simulation, but which nonetheless might be useful at some future point.

    /**
     * Calculates the natural logarithm of the Gamma function using the
     * Lanczos approximation with g=4.7421875 and n=15.
     *
     * @param n the input to the Gamma function.
     * @return the output of the Gamma function.
     */
    public static double lnGamma(double n){
        if (n <= 0) throw new IllegalArgumentException("n must be >0.");
        double[] cof={57.1562356658629235,-59.5979603554754912,
                14.1360979747417471,-0.491913816097620199,.339946499848118887e-4,
                .465236289270485756e-4,-.983744753048795646e-4,.158088703224912494e-3,
                -.210264441724104883e-3,.217439618115212643e-3,-.164318106536763890e-3,
                .844182239838527433e-4,-.261908384015814087e-4,.368991826595316234e-5};
        double x,tmp,y,ser;

        y=x=n;
        tmp = x+5.24218750000000000; // Rational 671/128
        tmp = (x+0.5)*log(tmp)-tmp;
        ser = 0.999999999999997092;
        for (int j=0; j < 14; j++) ser += cof[j]/++y;
        return tmp+log(2.5066282746310005*ser/x); //sqrt(2pi)
    }

    /**
     * Given distinct initial points ax and bx, searches in the downhill
     * direction (defined by the likelihood function as evaluated at the
     * initial points) and returns new points ax, bx, cx that bound the minimum
     * in a List.
     *
     * @param a lower initial bound.
     * @param b higher initial bound.
     * @param mu the mean number of mutations per generation.
     * @param M the number of samples.
     * @param etaTilda sample frequency spectrum.
     * @param aTildas switch controlling whether aTildas or aDeltas are used in
     *        likelihood calculations, with true corresponding to aTilda.
     * @return a List containing bounds to be inputted into a golden-section search.
     */
    public static List<Double> goldenBounds(final double a, final double b, double mu, int M, List<Integer> etaTilda, boolean aTildas) {
        final double GOLD=1.618034,GLIMIT=100.0,TINY=1.0e-20;
        double ax,bx,cx,fa,fb,fc;
        ax=a; bx=b;
        double fu;
        fa=-logLikelihood(mu, ax, M, etaTilda, aTildas);
        fb=-logLikelihood(mu, bx, M, etaTilda, aTildas);
        if (fb > fa) {
            double swap=ax; ax=bx; bx=swap;
            swap=fb; fb=fa; fa=swap;
        }
        cx=bx+GOLD*(bx-ax);
        fc=-logLikelihood(mu, cx, M, etaTilda, aTildas);
        while (fb > fc) {
            double r=(bx-ax)*(fb-fc);
            double q=(bx-cx)*(fb-fa);
            double u=bx-((bx-cx)*q-(bx-ax)*r)/
                    (2.0*SIGN(max(abs(q-r),TINY),q-r));
            double ulim=bx+GLIMIT*(cx-bx);
            if ((bx-u)*(u-cx) > 0.0) {
                fu=-logLikelihood(mu, u, M, etaTilda, aTildas);
                if (fu < fc) {
                    ax=bx;
                    bx=u;
                    List<Double> bounds = new ArrayList<>();
                    bounds.add(ax); bounds.add(bx); bounds.add(cx);
                    return bounds;
                } else if (fu > fb) {
                    cx=u;
                    List<Double> bounds = new ArrayList<>();
                    bounds.add(ax); bounds.add(bx); bounds.add(cx);
                    return bounds;
                }
                u=cx+GOLD*(cx-bx);
            } else if ((cx-u)*(u-ulim) > 0.0) {
                fu=-logLikelihood(mu, u, M, etaTilda, aTildas);
                if (fu < fc) {
                    double dum = u+GOLD*(u-cx);
                    bx=cx; cx=u; u=dum;
                }
            } else if ((u-ulim)*(ulim-cx) >= 0.0) {
                u=ulim;
            } else {
                u=cx+GOLD*(cx-bx);
            }
            ax=bx; bx=cx; cx=u;
        }
        List<Double> bounds = new ArrayList<>();
        bounds.add(ax); bounds.add(bx); bounds.add(cx);
        return bounds;
    }

    /**
     * An Euler-Maclaurin approximation which turns S(k=0, m, n, r) = S_0mn(r)
     * from a discrete sum into an integral with a closed form.
     *
     * @param r the basic reproduction number.
     * @param m the power to which r^-g is raised in the original sum.
     * @param n the power to which 1 - r^-g is raised in the original sum.
     * @return an approximated, but untruncated calculation of S_0mn(r).
     * @see <img src="./doc-files/S0mn.png"/>
     */
    public static double S_0mn(double r, double m, double n){
        BigInteger term1 = factorial(m-1).multiply(factorial(n));
        BigInteger term2 = factorial(m+n);
        BigInteger term = term1.divide(term2);
        double m1 = 0;
        double n1 = 0;
        if(m==0){m1=1;}
        if(n==0){n1=1;}
        return term.doubleValue()*log(r) + 0.5*(n1 + m1);
    }

    /**
     * Directly calculates n! using the BigInteger class to handle overflow.
     *
     * @param n the number whose factorial is to be calculated.
     * @return n!
     */
    public static BigInteger factorial(double n) {
        BigInteger result = BigInteger.valueOf(1);
        for (long factor = 2; factor <= n; factor++) {
            result = result.multiply(BigInteger.valueOf(factor));
        }
        return result;
    }

    /**
     * Approximates a factorial by using the Lanczos approximation of the Gamma
     * function.
     *
     * @param n the number whose factorial is to be approximated.
     * @return ~n!
     */
    public static double factorialG(double n){
        return exp(lnGamma(n+1));
    }
}
