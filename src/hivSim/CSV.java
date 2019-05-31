package hivSim;

import nr.ran.Poissondev;
import nr.ran.Ran;

import java.io.*;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

import static hivSim.MLE.*;
import static java.lang.StrictMath.*;

/**
 * This class handles input/output and integrates methods from all classes in
 * the program. main() will only call methods in this class.
 */
public class CSV {
    //TODO: Output random seeds
    //TODO: Profile program

    /**
     * @param i the number of significant figures desired.
     */
    public static void setTolerance(int i){
        MLE.setTolerance(i);
        MOM.setTolerance(i);
    }

    /**
     * Generates a specifiable number of phylogenies of VirusTrees and then
     * attempts to estimate the r0 used to generate the phylogeny by using
     * maximum likelihood estimation and the method of moments. Estimates and
     * standard deviations are then outputted into .CSV files. Raw data and
     * means are outputted in separate sheets, as are different sample sizes.
     *
     * @param N the number of Viruses at which to limit the phylogenetic
     *        simulation.
     * @param mu the mean number of mutations per generation.
     * @param GammaModel controls whether Viral reproduction follows the Gamma
     *        model
     * @param randL the initial random seed.
     * @param realizations the number of phylogenies to simulate for each r0.
     * @param tolerance the number of significant figures required.
     * @param aTildas switch controlling whether aTildas or cBars are used in
     *        likelihood calculations, with true corresponding to cTilda.
     */
    public static void simulate(int N, double mu, boolean GammaModel, long randL, int realizations, int tolerance, boolean aTildas) throws IOException {

        //initialize
        long runtime = System.currentTimeMillis();
        setTolerance(tolerance);
        long randL_old = randL; //store for output
        List<Double> rValues = rValues(1, 10, 0.1);
        List<Integer> MValues = MValues(8, 1024);
        List<List<List<Double>>> data = MValues.stream().mapToInt(x -> x).<List<List<Double>>>mapToObj(x -> new ArrayList<>()).collect(Collectors.toList()); //make data storage Lists for each M
        populateBinom(MValues.get(MValues.size() - 1)); //set binomial lookup table (array) to appropriate size

        //for (int M: MValues){
        //    FileWriter fw = new FileWriter(new File(M + "_" + mu + ".csv"));
        //    fw.close();
        //}

        //calculate
        for (double r0: rValues){
            //populateSkmn(2, 1024, 1024);
            for (int i = 1; i <= realizations; i++){

                System.out.printf("Seconds since last iteration: %f%n", (System.currentTimeMillis() - runtime)/1000.0);
                System.out.printf("r0: %f;  i: %d%n", r0, i);
                runtime = System.currentTimeMillis();

                Ran Random = new Ran(randL); randL = Random.int64(); Poissondev poisson = new Poissondev(r0, randL); //generate new random seed
                VirusTree Tree = new VirusTree(new Virus<>( null, 0));
                List<VirusTree> live = VirusTree.phylogenySim(Tree, N, r0, GammaModel, randL, poisson);
                for (int M: MValues){

                    //List<Double> etaTilda = toDouble(VirusTree.autoFoldedSFS(live, Tree, M, mu, poisson, randL, 0)); //ML training data
                    //etaTilda.add(r0);
                    //NN.appendData(etaTilda, M + "_" + mu);
                    System.out.printf("%d; ", M);
                    List<Double> calculate = calculate(M, r0, mu, randL, live, Tree, poisson, aTildas);
                    data.get(MValues.indexOf(M)).add(calculate); //add calculated values to data
                }
            }
        }

        //output
        for (int M: MValues){
            int index = MValues.indexOf(M);
            List<List<Double>> means = listCruncher(data.get(index), rValues);
            outputData(M, randL_old, realizations, 1, 10, 0.1, mu, data.get(index), "data" + M);
            outputData(M, randL_old, realizations, 1, 10, 0.1, mu, means,           "means" + M);
        }
    }

    public static List<Double> toDouble(List<Integer> intList){
        List<Double> listDouble = new ArrayList<>();
        for (Integer dx: intList){
            listDouble.add(dx.doubleValue());
        }
        return listDouble;
    }

    /**
     * Calculates estimates of r0 (and associated standard deviations in
     * log10 r0) given a List of live Viruses. Presently, this method uses two
     * versions of a maximum likelihood estimation (with and without an Euler-
     * McLaurin approximation) and the method of moments.
     *
     * Note that placeholder values are inserted after each set of values,
     * reserved for calculations of the sample standard deviation once all
     * realizations of the simulation have been handled.
     *
     * @param M number of samples to be taken from the live VirusTrees List.
     * @param r0 the basic reproduction number used in the simulation.
     * @param mu the mean number of mutations per generation.
     * @param randL the inherited random seed.
     * @param live the List of VirusTrees from which to sample.
     * @param Tree the founding VirusTree. //TODO: Multiple founders
     * @param poisson the inherited Poisson object used to generate a site-
     *        frequency spectrum.
     * @param aTildas switch controlling whether aTildas or cBars are used in
     *        likelihood calculations, with true corresponding to cTilda. 
     * @return a List of Doubles which contains, in order, r0, log10 r0, and
     *         data from the MLE, MLE_EM, and MOM. Each method's results
     *         include rHat, log10 rHat, log10 sHat, and a placeholder value
     *         for the sample standard deviation.
     */
    public static List<Double> calculate(int M, double r0, double mu, long randL, List<VirusTree> live, VirusTree Tree, Poissondev poisson, boolean aTildas){
        List<Integer> etaTilda = VirusTree.autoFoldedSFS(live, Tree, M, mu, poisson, randL, 0);
        double failures = etaTilda.get(etaTilda.size() - 1);
        etaTilda.remove(etaTilda.size() - 1);

        double log10e = log10(2.718281828459045235360287471); //conversion factor from sd(ln r) to sd(log10 r) is log10(e)
        List<Double> data = new ArrayList<>();

        double r_MLE = MLE.rHat(etaTilda, M, mu, aTildas);
        double s_hat_MLE = MLE.sigmaHat(r_MLE, M, mu, etaTilda, aTildas);

        double r_MLE_EM = MLE.rHatEM(M, mu, etaTilda);
        double s_hat_MLE_EM = MLE.sigmaHatEM(r_MLE_EM, M, etaTilda);

        double r_MOM = MOM.rHat(mu, M, etaTilda);
        double s_hat_MOM = MOM.sigmaHat(M, r_MOM, mu);

        data.add(r0);
        data.add(log10(r0));

        data.add(r_MLE);
        data.add(log10(r_MLE));
        data.add(log10e*s_hat_MLE);
        data.add(0.0);

        data.add(r_MLE_EM);
        data.add(log10(r_MLE_EM));
        data.add(log10e*s_hat_MLE_EM);
        data.add(0.0);

        data.add(r_MOM);
        data.add(log10(r_MOM));
        data.add(log10e*s_hat_MOM);
        data.add(0.0);

        data.add(failures);

        return data;
    }

    /**
     * Outputs all the elements of the given data set into a new .CSV file.
     * Parameters of the simulation will be included at the top of the file.
     *
     * @param M the number of samples taken.
     * @param randL the original random seed.
     * @param realizations the number of times a phylogeny was generated.
     * @param rMin the minimum r0 simulated.
     * @param rMax the maximum r0 simulated.
     * @param mu the number of mutations per generation.
     * @param data a List of Lists which stores data either raw calculate()
     *        data or averages of such data.
     * @param fileName the name of the output file.
     */
    public static void outputData(int M, long randL, int realizations, double rMin, double rMax, double r_delta, double mu, List<List<Double>> data, String fileName)throws FileNotFoundException{
        PrintWriter pw = new PrintWriter(new File(fileName + ".csv"));
        StringBuilder sb = new StringBuilder();

        sb.append("sample size,").append(M).append('\n');
        sb.append("realizations,").append(realizations).append('\n');
        sb.append("r0-,").append(rMin).append('\n');
        sb.append("r0+,").append(rMax).append('\n');
        sb.append("r0_d,").append(r_delta).append('\n');
        sb.append("mu,").append(mu).append('\n');
        sb.append("seed,").append(randL).append('\n');
        sb.append('\n');
        sb.append('\n');
        sb.append('\n');

        sb.append("r0,log10 r0,MLE r,MLE log10 r,MLE s^,MLE s,MLE_EM r,MLE_EM log10 r,MLE_EM s^,MLE_EM s,MOM r,MOM log10 r,MOM s^,MOM s,Failures,Successes,Failure Rate" + '\n');
        int size = data.get(0).size(); //# parameters to calculate
        for (List<Double> x: data){
            for (int i = 0; i < size-1; i++){
                sb.append(x.get(i).toString()); sb.append(',');
            }
            sb.append(x.get(size-1).toString()); sb.append('\n'); //need to press enter instead of , after last element
        }
        pw.write(sb.toString());
        pw.close();
    }

    /**
     * Turns raw data generated by calculate() into a means and standard
     * deviations.
     *
     * @param data a List of Lists which stores data generated by calculate().
     * @param rValues the inherited List of r0 values which were simulated.
     * @return a List of Lists whose outer index is the index of a given r0 in
     *         rValues and whose inner index corresponds to the mean or
     *         standard deviation of a variable from calculate().
     */
    public static List<List<Double>> listCruncher(List<List<Double>> data, List<Double> rValues){
        List<List<Double>> means = vectorMeans(data, rValues);
        int i = 0;
        for (List<Double> dt: means){
            dt.set(5, stdDev(data, rValues.get(i), 3)); //TODO: generalize these inputs
            dt.set(9, stdDev(data, rValues.get(i), 7));
            dt.set(13, stdDev(data, rValues.get(i), 11));
            dt.add(15, (double) (data.size() / rValues.size())); //number of successes
            dt.add(16, dt.get(14) / (dt.get(14) + dt.get(15))); //failure rate calculation
            i++;
        }
        return means;
    }

    /**
     * Calculates means for all variables stored in data (as produced by
     * calculate()) for each given r0.
     *
     * @param data a List of Lists which stores data generated by calculate().
     * @param rValues an inherited List of r0 values by which to index each set
     *        of averages.
     * @return a List of Lists, where each outer index corresponds to a
     *         simulation r0 and each inner index corresponds to the mean of a
     *         variable over all realizations of the simulation r0.
     */
    public static List<List<Double>> vectorMeans(List<List<Double>> data, List<Double> rValues){
        List<List<Double>> means = new ArrayList<>();
        for (double r0: rValues){
            means.add(listAverage(data, r0));
        }
        return means;
    }

    /**
     * Calculate averages for all variables stored in data (as produced by
     * calculate()) for a single given r0.
     *
     * @param data a List of Lists which stores data generated by calculate().
     * @param r0 the basic reproduction number.
     * @return a List containing averages for all values matching a given
     *         simulation r0. Note that the sample standard deviation average
     *         will return from here as 0, but will be calculated in the next
     *         method call.
     */
    public static List<Double> listAverage(List<List<Double>> data, double r0){
        List<Double> means = new ArrayList<>();
        for (int x = 0; x < data.get(0).size() - 1; x++){ //-1 because we sum the last column instead of averaging
            means.add(indexAverage(data, r0, x));
        }
        means.add(indexSum(data, r0, 14));
        return means;
    }

    /**
     * Calculates the sum of a variable from the data output as generated
     * by calculate().
     *
     * @param data a List of Lists which stores data generated by calculate().
     * @param r0 the basic reproduction number.
     * @param index the (inner) index.
     * @return the sum of a variable stored at the given inner index in the
     *         data vector for a given r0.
     */
    public static double indexSum(List<List<Double>> data, double r0, int index){
        List<List<Double>> temp = new ArrayList<>();
        for (List<Double> x: data){ //look at all realizations
            if (x.get(0) == r0){temp.add(x);} //extract only the ones corresponding to r0
        }
        double sum = 0;
        for (List<Double> dx: temp){
            sum = sum + dx.get(index);
        }
        return sum;
    }

    /**
     * Calculates the average of a variable from the data output as generated
     * by calculate().
     *
     * @param data a List of Lists which stores data generated by calculate().
     * @param r0 the basic reproduction number.
     * @param index the (inner) index.
     * @return the average of a variable stored at the given inner index in the
     *         data vector for a given r0.
     */
    public static double indexAverage(List<List<Double>> data, double r0, int index){
        List<List<Double>> temp = new ArrayList<>();
        for (List<Double> x: data){
            if (x.get(0) == r0){temp.add(x);}
        }
        double sum = 0;
        for (List<Double> dx: temp){
            sum = sum + dx.get(index);
        }
        return sum/temp.size();
    }

    /**
     * Calculates the sample standard deviation (using n-1 instead of n) in r0
     * for a variable in the data output as generated by calculate().
     *
     * @param data a List of Lists which stores data generated by calculate().
     * @param r0 the basic reproduction number.
     * @param index the (inner) index.
     * @return the sample standard deviation of a variable stored at the given
     *         inner index in the data vector for a given r0.
     */
    public static double stdDev(List<List<Double>> data, double r0, int index){
        double xbar = indexAverage(data, r0, index);
        List<List<Double>> temp = new ArrayList<>();
        for (List<Double> x: data){
            if (x.get(0) == r0){temp.add(x);}
        }
        double sum = 0;
        for (List<Double> dx: temp){
            sum = sum + pow((xbar - dx.get(index)), 2);
        }
        sum = sum/(temp.size() - 1);
        return pow(sum, 0.5);
    }


    /**
     * Generates a monotonically ordered List of numbers from min to max,
     * exclusive, where index k stores the value (min + delta)^k+1.
     *
     * @param min defines the smallest element of the List: rmin + delta. Must
     *        be greater than 0.
     * @param max the (unreached) cut-off value. Must be larger than min.
     * @param delta the logarithmic difference between each term.
     * @return an ordered List from rmin + delta to the largest power of
     *         min + delta smaller than max, where index k stores the
     *         value (min + delta)^k+1.
     */
    public static List<Double> rValues(double min, double max, double delta){
        if (min < 0){throw new IllegalArgumentException("Minimum value must be >0.");}
        if (max < min){throw new IllegalArgumentException("Maximum value must be larger than minimum value.");}
        List<Double> rValues = new ArrayList<>();
        double k = 1;
        while (pow((min + delta), k) < max){
            rValues.add(pow((min + delta), k));
            k++;
        }
        return rValues;
    }

    /**
     * Generates a List storing consecutive powers of 2.
     *
     * @param min the smallest power of 2 to include.
     * @param max the largest power of 2 to include.
     * @return an ordered List containing consecutive powers of 2.
     */
    public static List<Integer> MValues(int min, int max){
        List<Integer> MValues = new ArrayList<>();
        for (int n=min; n<=max; n=2*n) MValues.add(n);
        return MValues;
    }
}
