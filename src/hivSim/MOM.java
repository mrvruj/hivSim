package hivSim;

import java.util.ArrayList;
import java.util.List;

import static hivSim.CSV.rValues;
import static java.lang.Math.abs;
import static java.lang.Math.floor;
import static java.lang.StrictMath.pow;

/**
 * Contains methods related to calculating the maximum likelihood estimate and
 * standard deviations of the basic reproduction number.
 *
 * Contains a global class variable 'tolerance' which controls the number of
 * significant figures desired; the default value is 7. The method
 * CSV.setTolerance() which calls MOM.setTolerance() and MLE.setTolerance(),
 * should be used to set this value.
 *
 * The constant value used for the arbitrary coefficients a_m used throughout
 * calculations in this class.
 */
public class MOM {
    private static int tolerance;
    private static double tol = pow(10, -tolerance);

    /**
     * @param i the number of significant figures desired.
     */
    static void setTolerance(int i){
        tolerance=i;
    }

    /**
     * Calculates an estimate for r0 using the method of moments, evaluating
     * only the r0's in rValues.
     *
     * @param mu the expected number of mutations per generation.
     * @param M the number of samples taken.
     * @param etaTilda an ordered List of values representing the folded site-
     *        frequency spectrum generated from the phylogenetic
     *        simulation.
     * @return an estimate of r0 according to the method of moments.
     * @see <img src="./doc-files/rHatMOM.png"/>
     */
    public static double rHat(double mu, int M, List<Integer> etaTilda){
        double bound = gridSearch(rValues(1, 10, 0.1), mu, M, etaTilda);
        return goldenSection(1.001, bound, 15, M, mu, etaTilda);
    }

    /**
     * Calculates a double value representing an estimate of the standard
     * deviation in ln rHat for the method of moments at a particular value
     * of r0.
     *
     * @param M the number of samples taken.
     * @param rHat the predicted r0 according to the method of moments.
     * @param mu the expected number of mutations per generation.
     * @return an estimate of the standard deviation for the method of moments
     *         at a particular value of r0.
     * @see <img src="./doc-files/sdMOM.png"/>
     */
    public static double sigmaHat(int M, double rHat, double mu){
        double h = h(M, rHat, mu);
        double denom = rHat*hPrime(M, mu, rHat);
        double expression = Math.sqrt(h);
        return abs(expression/denom);
    }

    /**
     * Uses a golden-section search to find the method of moments estimate of
     * the basic reproduction number.
     *
     * @param ax initial lower bound.
     * @param bx initial midpoint.
     * @param cx initial upper bound.
     * @param M sample size.
     * @param mu mean number of mutations per generation.
     * @param etaTilda folded site frequency spectrum.
     * @return the method of moments for the basic reproduction number as
     *         estimated by golden section minimization.
     * @see <img src="./doc-files/goldenSection.png"/>
     */

    private static double goldenSection(double ax, double bx, double cx, int M, double mu, List<Integer> etaTilda)  {
        final double R=0.61803399,C=1.0-R;
        double xmin;

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
        double f1= abs(etaSum(M, etaTilda) - h(M, x1,mu));
        double f2= abs(etaSum(M, etaTilda) - h(M, x2,mu));
        while (abs(x3-x0) > tol*(abs(x1)+abs(x2))) {
            if (f2 < f1) {
                double dum = R*x2+C*x3;
                x0=x1; x1=x2; x2=dum;
                f1=f2;
                f2 = abs(etaSum(M, etaTilda) - h(M,x2,mu));
            } else {
                double dum = R*x1+C*x0;
                x3=x2; x2=x1; x1=dum;
                f2=f1;
                f1= abs(etaSum(M, etaTilda) - h(M,x1,mu));
            }
        }
        if (f1 < f2){xmin=x1;}
        else        {xmin=x2;}
        return xmin;
    }

    /**
     * Calculates and sums h_m(r) from m = 1 to m = M-2. h_m(r) is the product
     * of a_m, mu, and C_bar(r) and represents the expectation of a similar sum
     * of eta. An array of these sums for different values of r will then be
     * compared to the data generated from the phylogenetic simulation.
     *
     * @param M the number of samples taken.
     * @param r0 the basic reproduction number for which to generate the sum.
     * @param mu the expected number of mutations per generation.
     * @return a double representing the sum of a_m * h_m for m = 2 to m = M-2.
     * @see <img src="./doc-files/h(r).png"/>
     */
    public static double h(int M, double r0, double mu){
        double term;
        double sum = 0;
        for (int m = 2; m <= M-2; m++){
            term = mu*MLE.aDelta(r0, M,m);
            sum += term;
        }
        return sum;
    }

    /**
     * Calculates and sums the derivative of h_m(r) from m = 1 to m = M-2.
     * This method is used to calculate the estimated standard deviation for
     * the method of moments.
     *
     * @param M the number of samples taken.
     * @param mu the expected number of mutations per generation.
     * @param r0 the basic reproduction number for which to generate the sum.
     * @return a double representing the derivative of the sum of a_m * h_m
     *         for m = 2 to m = M-2.
     * @see <img src="./doc-files/hPrime.png"/>
     */
    static double hPrime(int M, double mu, double r0){
        double sum = 0;
        for (int m = 2; m <= M-2; m++){
            double term = mu*MLE.aDeltaPrime(r0, M, m);
            sum += term;
        }
        return sum;
    }

    /**
     * Sums the elements in a List of etaTilda's from m = 2 to m = M_tilda - 1.
     *
     * @param M the number of samples taken.
     * @param etaTilda an ordered List of values representing the folded site-
     *        frequency spectrum generated from the phylogenetic
     *        simulation.
     * @return a double representing the sum of etaTilda_m for m = 2 to
     *         m = M-2.
     * @see <img src="./doc-files/etaSum.png"/>
     */
    static double etaSum(int M, List<Integer> etaTilda){
        double sum=0;
        for (int m = 2; m <= floor(M/2.0); m++){
            sum += etaTilda.get(m-1); //indices start at 0
        }
        return sum;
    }

    /**
     * Uses a grid to identify the method of moments estimate of the basic
     * reproduction number.
     *
     * @param rValues an ordered List of Doubles which represents the values of
     *        r0 to be tested.
     * @param mu the expected number of mutations per generation.
     * @param M the number of samples taken.
     * @param etaTilda an ordered List of values representing the folded site-
     *        frequency spectrum generated from the phylogenetic
     *        simulation.
     * @return the value in the input grid which produces the smallest
     *         abs(h(r0) - etaTilda) as a function of r0.
     */
    static double gridSearch(List<Double> rValues, double mu, int M, List<Integer> etaTilda){
        List<Double> grid = new ArrayList<>();
        for (double r0: rValues){
            grid.add((abs(h(M, r0, mu) - etaSum(M, etaTilda))));
        }
        int smallest = 0;
        for (int i=1; i<grid.size(); i++){
            if (grid.get(i) < grid.get(smallest)){smallest = i;}
        }
        return rValues.get(smallest);
    }
}
