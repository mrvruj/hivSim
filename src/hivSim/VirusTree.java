package hivSim;

import nr.ran.*;
import java.util.*;
import java.util.stream.Collectors;
import static java.lang.Math.floor;
import static java.lang.Math.max;
import static java.lang.StrictMath.pow;

//TODO: Add a method remove()
//TODO: replace etaTilda List with a double[] for cleanliness, i.e. fixing functions using m-1 instead of m

/**
 * An instance represents a tree of Viruses. Each Virus in the VirusTree
 * may have only one parent Virus. Each Virus in the VirusTree is the child
 * of the Virus it is immediately descended from. A VirusTree may or may not
 * have a parent VirusTree.
 *
 * The founder of the VirusTree does not have a parent Virus. Each element of
 * children is descended from the founder. Children is non-null but will be an
 * empty List if the founder is a leaf.
 *
 * This class contains only methods relating to the generation of a VirusTree
 * phylogeny.
 **/
public class VirusTree {
    private Virus founder;
    private List<VirusTree> children;

    /**
     * Constructor: a new VirusTree with founder p and no children.
     *
     * @param p the founder of this VirusTree.
     */
    public VirusTree(Virus p){
        founder = p;
        children = new ArrayList<>();
    }

    /**
     * Constructor: a new VirusTree that is a copy of tree p.
     * Tree p and its copy have no VirusTrees in common.
     *
     * @param p the VirusTree to be copied.
     */
    public VirusTree(VirusTree p){
        founder = p.founder;
        children = new ArrayList<>();
        for (VirusTree dt: children){
            children.add(new VirusTree(dt));
        }
    }

    /** @return the founder of this VirusTree. */
    public Virus getFounder() {
        return founder;
    }

    /** @return A List of children VirusTree objects. */
    public List<VirusTree> getChildren() {
        return this.children;
    }

    /** @return the data stored in the (Virus) founder of this VirusTree. */
    @Override
    public String toString() {
        return this.getFounder().toString();
    }

    /**
     * Searches this VirusTree recursively for a VirusTree whose founder is p.
     *
     * @param p the (Virus) founder of a VirusTree.
     * @return the VirusTree object in this VirusTree whose founder is p or
     *         null if p is not in the VirusTree.
     *
     * TODO: Multiple Founders
     */
    public VirusTree getTree(Virus p) {
        if (founder == p) return this;  //base case
        for (VirusTree dt: children) { //recursive case
            VirusTree search= dt.getTree(p);
            if (search != null) return search;
        }
        return null; // p is not in the tree
    }

    /**
     * Recursively checks if this VirusTree contains p iff the founder of this
     * VirusTree is p or if one of the founder's descendants contains p.
     *
     * @param p the Virus to be searched for.
     * @return True iff this VirusTree contains p.
     */
    public boolean contains(Virus p) {
        if (founder == p) return true;  //base case
        for (VirusTree dt: children) { //recursive case
            if (dt.contains(p)) return true;
        }
        return false;
    }

    /**
     * Adds Virus c to this VirusTree as a child of Virus p.
     *
     * @param p the parent virus. Must be non-null.
     * @param c the child virus. Must be non-null.
     * @return the VirusTree whose founder is the new child.
     */
    public VirusTree add(Virus p, Virus c) {
        if (p == null || c == null){throw new IllegalArgumentException("Parent and child Viruses must be non-null.");}
        VirusTree parent = getTree(p);
        VirusTree child = new VirusTree(c);
        parent.children.add(child);
        return child;
    }

    /**
     * Recursively counts the number of Viruses in this VirusTree.
     *
     * @return the number of Viruses in this VirusTree.
     */
    public int size() {
        int size = 1;
        if(this.getChildren().size() == 0){ //base case
            return 1;
        }
        else{
            for (VirusTree dt: children){  //recursive case
                size = size + dt.size();
            }
        }
        return size;
    }

    /**
     * Removes all children from VirusTree v.
     *
     * @param v the VirusTree whose children will be removed.
     */
    public static void removeChildren(VirusTree v){
        v.children.removeAll(v.children);
    }

    /**
     * Recursively add all VirusTrees descended from this VirusTree without
     * any children to a List.
     *
     * @return all childless VirusTrees (i.e. leafs on the tree) descended
     *         from this VirusTree.
     */
    public List<VirusTree> getChildless(){
        List<VirusTree> list = new LinkedList<>();
        if (this.getChildren().size() == 0){ //base case
            list.add(this);
        }
        for (VirusTree dt: children){   //recursive case
            list.addAll(dt.getChildless());
        }
        return list;
    }

    /**
     * Adds a random (according to a Poisson distribution) number of daughters
     * to parent in the given VirusTree with a random(according to a gamma
     * distribution) burstday.
     *
     * @param parent the Virus whose VirusTree to which children are added to.
     *        Must be non-null.
     * @param r0 the basic reproduction number. Must be >0.0.
     * @param Gamma true if using non-deterministic birth/burst times,
     *        according to the Gamma model.
     * @param gammaN the shape parameter of the gamma distribution.
     * @param gammaL the rate parameter of the gamma distribution. TODO: use gammaParamters()
     * @param randomSeed the random number inherited from the simulation.
     * @param poisson the Poisson object inherited from the simulation.
     */
    public void addRandom(Virus parent, double r0, boolean Gamma, double gammaN, double gammaL, Poissondev poisson, long randomSeed){
        if (r0 < 0){
            throw new IllegalArgumentException("r0 must be >0.");
        }
        if (parent == null){
            throw new IllegalArgumentException("Virus parent must be non-null.");
        }
        int numD = poisson.dev(r0); // number of offspring
        for (double x=0; x< numD; x++){ // add daughters
            double rand2 = 2; // fixed burstday if delta
            if (Gamma){
                Gammadev gamma = new Gammadev(gammaN, gammaL, randomSeed);
                rand2 = gamma.dev(); //random burstday if gamma
            }
            Virus child = new Virus(parent.getBurstday() + rand2, parent.getBurstday() + rand2); //add randomly generated number to burstday of parent
            this.add(parent, child);
        }
    }

    /**
     * Randomly (according to a gamma and poisson distribution) adds generations
     * of children to a given founding VirusTree. Stores the final generation
     * in a List. If the lineage of the founder goes extinct, flush the List
     * and restart the simulation until the size of the List is maxVirus.
     *
     * @param foundingVirusTree the founding VirusTree. Must be non-null.
     * @param maxVirus the number of (live) Viruses at which the
     *        simulation terminates. Must be >0.0.
     * @param r0 the basic reproduction number. Must be >0.0.
     * @param GammaModel a boolean value, where true means the simulation
     *        is using the gamma model to generate offspring burstdays.
     * @param randL the inherited random seed.
     * @param poisson the inherited poisson number generator.
     * @return a list of VirusTrees whose founders are the last Viruses
     *         concurrently live in the phylogeny.
     *
     * TODO: Multiple Founders - convert foundingVirusTree -> List<VirusTree>
     *  and use for-each to addRandom(). Gonna have to pay attention to the
     *  sort. Maybe use if-statements?
     *  (if tree1.contains(live.get(0)){tree1.addRandom()}
     */
    public static List<VirusTree> phylogenySim(VirusTree foundingVirusTree, int maxVirus, double r0, boolean GammaModel, long randL, Poissondev poisson){
        removeChildren(foundingVirusTree);
        if (r0 < 0){
            throw new IllegalArgumentException("r0 must be >0.");
        }
        if (maxVirus < 0){
            throw new IllegalArgumentException("maxVirus must be >0.");
        }
        try{
            List<VirusTree> live = foundingVirusTree.getChildless();
            while (live.size() < maxVirus){
                if (live.isEmpty()){ // handle extinct virus
                    throw new Exception("Virus went extinct.");
                }
                foundingVirusTree.addRandom(live.get(0).getFounder(), r0, GammaModel, 69.4, 34.7, poisson, randL); //add daughters to the oldest live virus
                for (VirusTree dt : live.get(0).getChildren()){ //add newly born children to live
                    merge(live, dt); //place them in live according to burstday
                }
                live.remove(0); //remove burst Virus from live
            }
            return live;
        }
        catch(Exception e){
            Ran uniform = new Ran(randL); //generate new random seed
            return phylogenySim(foundingVirusTree, maxVirus, r0, GammaModel, uniform.int64(), poisson);
        }
    }

    /**
     * Merges the VirusTree child into the given list. Starts the search from
     * the end because most newly born VirusTrees will be the youngest in the
     * list.
     *
     * @param list the list of VirusTrees into which a new element is to be
     *        merged.
     * @param child the VirusTree to be merged.
     */
    public static List<VirusTree> merge(List<VirusTree> list, VirusTree child){
        int i = list.size() - 1;
        while (i >= 0){
            if (list.get(i).isOlderThan(child)){
                list.add(i, child);
                return list; //return to end method call when child is added
            }
            i--;
        }
        return list;
    }

    /**
     * Checks if this VirusTree has a smaller (corresponding to older) burstday
     * than tree's burstday.
     *
     * @param tree the VirusTree to be compared to this VirusTree.
     * @return true iff this VirusTree is older than tree.
     */
    public boolean isOlderThan(VirusTree tree){
        return this.getFounder().getBurstday() <= tree.getFounder().getBurstday();
    }

    /**
     * Generates parameters of the Gamma distribution given a mean and
     * standard deviation. Note the alpha, beta parameterization is used
     * instead of the k, theta. The mean is alpha/beta and the variance is
     * alpha/beta^2.
     *
     * @param mean the mean of the distribution.
     * @param std the standard deviation of the distribution.
     * @return a List where the first element is alpha and the second is beta.
     * //TODO: add @see
     */
    public static List<Double> gammaParameters(double mean, double std){
        List<Double> parameters = new ArrayList<>();
        double var = pow(std, 2);
        double mu2 = pow(mean, 2);
        double betaRadical = mean/var;
        double alphaRadical = var*mu2;
        parameters.add(pow(alphaRadical, 1.0/3));
        parameters.add(pow(betaRadical, 1.0/3));
        return parameters;
    }

    /**
     * Randomly (according to a uniform distribution) sample, without
     * replacement, n Viruses from a List of VirusTrees. Note that the random
     * number generator in this method uses Random() from java.util rather
     * than the corresponding method in Numerical Recipes because of
     * compatibility issues with Collections.shuffle(List<T>, Random x).
     *
     * @param live the List of VirusTrees to randomly sample from.
     * @param M the number of Viruses to be sampled.
     * @param randL inherited random seed.
     * @return a List containing n random Viruses from List live.
     */
    public static List<Virus> sample(List<VirusTree> live, double M, long randL) {
        List<Virus> sample = new LinkedList<>();
        Collections.shuffle(live, new Random(randL)); //TODO: Find way to inherit RNG from NR
        for (int i = 0; i < M; ++i){
            sample.add(0, live.get(i).getFounder());
        }
        return sample;
    }

    /**
     * Recursively finds and returns the ancestry of Virus c (including c) in a
     * List of Viruses.
     *
     * @param c the Virus whose ancestry we want.
     * @return the route the Virus took to get from the founder of this VirusTree
     *         to child c or null if c is not in the VirusTree.
     */
    public List<Virus> routeTo(Virus c) {
        if (c == founder){  //base case
            LinkedList<Virus> list = new LinkedList<>();
            list.addFirst(c);
            return list;
        }
        else {
            for (VirusTree dt: children){  //recursive case
                List<Virus> list = dt.routeTo(c);
                if (list != null){
                    list.add(0, founder);
                    return list;
                }
            }
        }
        return null;
    }

    /**
     * Given a sample of n Viruses, return a Map with each Virus in a VirusTree
     * assigned a Long representing how many Viruses in the sample are
     * descended from it. First finds the ancestry of all Viruses in the sample
     * and adds all found Viruses to a master List. Each occurrence of a given
     * Virus in the List raises it's running count (stored in the Map) by 1.
     *
     * TODO: Multiple Founders
     *
     * @param sample a List of Viruses from a common VirusTree.
     * @return a Map from all the Viruses in a VirusTree to Longs representing
     *         the number of direct descendant each Virus has in the sample.
     */
    public static Map sampleAncestryCounter(List<Virus> sample, VirusTree founder){
        List<Virus> ancestryHolder = new LinkedList<>();
        for (Virus dt: sample){
            ancestryHolder.addAll(founder.routeTo(dt));
        }
        return ancestryHolder.stream().collect(Collectors.groupingBy(e -> e, Collectors.counting()));
    }

    /**
     * Generates and returns the ancestor frequency spectrum as an array, where
     * element [i] is the count of ancestors with i+1 descendants in the sample.
     *
     * This method first counts the number of direct descendants each Virus in
     * the founding VirusTree has in the given sample and stores this information
     * in a Map<Virus, Long>. This Map is converted to a Map<Long, Long> taking
     * the Longs from the previous Map and attaching them to Longs representing
     * how many times they occur in the previous Map. This second Map is then
     * converted into an array, where element [i] is the count of ancestors with
     * i+1 descendants in the sample.
     *
     * TODO: Multiple Founders
     *
     * @param sample a List of Viruses from a common VirusTree.
     * @param founder is the founding VirusTree.
     * @return a List where the element [i] represents the number of ancestors
     *         with i+1 descendants in the sample.
     */
    public static List<Integer> ancestorFrequencySpectrum(List<Virus> sample, VirusTree founder) {
        int sampleSize = sample.size();
        Map<Virus, Long> ancestry = sampleAncestryCounter(sample, founder); //count how many times each Virus appears
        Map<Long, Long> AFS_map = ancestry.values().stream().collect(Collectors.groupingBy(e -> e, Collectors.counting()));

        List<Integer> AFS = new ArrayList<>(sampleSize);
        for (long i = 1; i <= sampleSize; ++i){ //turn above Maps into an array
            if (AFS_map.get(i) != null){AFS.add(AFS_map.get(i).intValue());}
            else{AFS.add((0));}
        }
        return AFS;
    }

    /**
     * Generates a site frequency spectrum by randomly (according to
     * a Poisson distribution and an empirical mutation rate) generating
     * mutations about a given ancestor frequency spectrum. In practice, a
     * site frequency spectrum is not observable because the founding
     * sequence is usually unknown.
     *
     * @param A a List where element [i] represents the number of ancestors
     *        with i+1 descendants in the sample.
     * @param mu the mutation rate per generation.
     * @param poisson a random Poisson number generator inherited from the
     *        simulation.
     * @return a List form of the site frequency spectrum where element [i]
     *         represents the count of i+1 letters in any Virus in the
     *         sample differing from the corresponding founder letter. Element
     *         [i] is a Poisson-distributed random variable centered on mu
     *         times the ancestor frequency spectrum.
     */

    public static List<Integer> siteFrequencySpectrum(List<Integer> A, double mu, Poissondev poisson){
        int M = A.size();
        List<Integer> eta = new ArrayList<>(M);
        for (int m = 0; m < M; m++){
            eta.add(m, poisson.dev(mu*A.get(m)));
        }
        return eta;
    }

    /**
     * Folds a given site frequency spectrum (eta) to make it agnostic to
     * whether a letter deviating from the majority letter corresponds to a
     * mutation or non-mutation. Generated by combining elements corresponding
     * to m or M-m minority letters, where M is the sample size.
     *
     * @param eta a List representing the site frequency spectrum, where
     *        element [i] represents the count of i+1 letters in any Virus in
     *        the sample differing from the corresponding founder letter.
     * @return a List form of the folded site frequency spectrum, where
     *         element [i] is the count of alignment columns where the
     *         number of minority letters is i+1.
     * @see <img src="./doc-files/foldedSiteFrequencySpectrum.png"/>
     */
    public static List<Integer> foldedSiteFrequencySpectrum(List<Integer> eta){
        double M = eta.size();
        int mTilda = (int) floor((M)/2.0);
        List<Integer> etaTilda = new ArrayList<>();
        for (int m = 0; m < mTilda - 1 ; m++){ //handle all indices except final index
            etaTilda.add(eta.get(m) + eta.get((int) M - 2 - m));
        }
        if (M % 2 == 0){ //handle final index if M is even
            etaTilda.add(eta.get(mTilda - 1));
        }
        else{ //handle final index if M is odd
            etaTilda.add(eta.get(mTilda) + eta.get((int) M - 2 - mTilda));
        }
        return etaTilda;
    }

    /**
     * Generates a folded site-frequency spectrum from a List of live
     * VirusTrees. If the sum of the elements of the SFS from 1 to floor(M/2)
     * is 0, resamples the List and generates a new SFS until a productive
     * SFS is found. WARNING: Last index counts number of failed attempts.
     *
     * @param live the List of VirusTrees to randomly sample from.
     * @param founder = the founding VirusTree. TODO: Multiple founders
     * @param M the number of Viruses to be sampled.
     * @param randL inherited random seed.
     * @param mu the mean number of mutations per generation.
     * @param poisson the inherited Poisson distribution object.
     * @param randL inherited random number.
     * @param failures counts the number of times this function goes through
     *        it's try-catch block. Should always be 0 outside of the single
     *        internal recursive call in the catch block.
     * @return a productive folded site-frequency spectrum with the last index
     *         counting the number of failed attempts preceded this productive
     *         SFS.
     */
    public static List<Integer> autoFoldedSFS(List<VirusTree> live, VirusTree founder, int M, double mu, Poissondev poisson, long randL, int failures) {
        try {
            List<Virus> sample = sample(live, M, randL);
            List<Integer> AFS = ancestorFrequencySpectrum(sample, founder);
            List<Integer> Eta = siteFrequencySpectrum(AFS, mu, poisson);
            List<Integer> etaTilda = foldedSiteFrequencySpectrum(Eta);
            double sum = 0;
            for (int i = 1; i < floor(M / 2.0); i++) {
                sum = sum + etaTilda.get(i);
            }
            if (sum == 0) {
                failures++;
                throw new Exception("No usable information in this SFS.");
            }
            etaTilda.add(failures);
            return etaTilda;
        } catch (Exception e) {
            Ran uniform = new Ran(randL); //uniform distribution
            return autoFoldedSFS(live, founder, M, mu, poisson, uniform.int64(), failures);
        }
    }

    //Below this point are functions which are not used in the main body of the
    //simulation, but which nonetheless might be useful at some future point.

    /**
     * Recursively add all Viruses descended from the founder of VirusTree p
     * to a List. The founder of p is the first element of this list, but the
     * list is not otherwise ordered.
     *
     * @return a List containing all Viruses in VirusTree p, including the
     *         founder of p.
     */
    public static List<Virus> getAll(VirusTree p){
        List<Virus> list = new LinkedList<>();
        list.add(p.founder);                    //base case
        for (VirusTree dt: p.getChildren()){    //recursive case
            list.addAll(getAll(dt));
        }
        return list;
    }

    /**
     * Finds the youngest common ancestor of Viruses child1 and child2.
     *
     * @param child1 one of two Viruses whose common ancestor is required.
     * @param child2 one of two Viruses whose common ancestor is required.
     * @return the founder Virus of the smallest subtree of this VirusTree
     *         that contains child1 and child2. If either child is null
     *         or is not in this VirusTree, return null.
     */
    public Virus sharedAncestorOf(Virus child1, Virus child2) {
        if (child1 == null || child2 == null){return null;}
        List<Virus> route1 = routeTo(child1);
        List<Virus> route2 = routeTo(child2);
        for (int g = 0; g <= route1.lastIndexOf(child1); g++){
            if (getTree(child1).contains(child2)){return child1;} //direct ancestry
            if (getTree(child2).contains(child1)){return child2;}
            if (route2.get(g) != route2.get(g)){return route1.get(g-1);}
        }
        return null;
    }

    /**
     * Searches this VirusTree recursively for the Virus whose child is c.
     *
     * @param c the Virus whose parent is needed.
     * @return the immediate (Virus) parent of c or null if c is the founder
     *         or not in this VirusTree.
     */
    public Virus getParent(Virus c) {
        for (VirusTree dt: children) {     //base case
            if (dt.founder == c) return founder;
        }
        for (VirusTree dt: children) {     //recursive case
            Virus parent= dt.getParent(c);
            if (parent != null) return parent;
        }
        return null; //c is the founder or not in the tree
    }

    /**
     * Recursively calculates the inclusive number of nodes between the founder
     * to the Virus p. Note: genOf(founder) is 0. If p is a child of this VirusTree,
     * then genOf(p) is 1.
     *
     * @param p the Virus whose generation is required.
     * @return the depth at which p occurs in this VirusTree, or -1
     *         if p is not in the VirusTree.
     */
    public int genOf(Virus p) {
        int gen;
        if (getTree(p) == null){return -1;}
        if (getParent(p) == null){  //base case
            return 0;
        }
        else {                      //recursive case
            gen = 1 + genOf(getParent(p));
        }
        return gen;
    }

    /**
     * Recursively finds all Viruses descended from this VirusTree in
     * generation g and add them to a List. Note g = 0 returns the founder
     * of this VirusTree.
     *
     * @param g the generation to be returned. Must be >= 0.
     * @return a List containing entire width in generation g.
     */
    public List<Virus> getGen(int g) {
        if (g < 0){
            throw new IllegalArgumentException("Generation number must be >=0.");
        }
        List<Virus> list = new LinkedList<>();
        if (this.genOf(founder) == g){  //base case
            list.add(this.founder);
        }
        for (VirusTree dt: children){   //recursive case
            list.addAll(dt.getGen(g-1));
        }
        return list;
    }

    /**
     * Recursively calculates the number of VirusTrees at a given generation of
     * this VirusTree's descendants.
     *
     * @param g the generation of this VirusTree to be measured. Must be >=0.
     * @return the width of this VirusTree at generation g, i.e. the number of
     *         VirusTrees that occur at generation g, where the generation of
     *         the founder is 0.
     */
    public int widthAtGen(int g) {
        if (g < 0){
            throw new IllegalArgumentException("Generation number must be >=0.");
        }
        int width = 0;
        if (g == 0){    //base case
            return 1;
        }
        else{     //recursive case
            for (VirusTree dt: children){
                width = width + dt.widthAtGen(g-1);
            }
        }
        return width;
    }

    /**
     * Recursively searches for the maximum depth, i.e. the longest path from
     * the founder to a leaf, of this VirusTree.
     *
     * @return the maximum depth of this VirusTree. If this VirusTree is a leaf,
     *         return 0.
     */
    public int maxGen() {
        int maxGen= 0;
        for (VirusTree dt: children) {
            maxGen = max(maxGen, dt.maxGen() + 1);
        }
        return maxGen;
    }

    /**
     * Calculates the generation of the last descendant of Virus c and compares
     * that value with the generation of the last descendant of Virus c's
     * oldest ancestor.
     *
     * @param c the Virus to be examined.
     * @return True if Virus c doesn't have any children in it's oldest
     *         ancestor's youngest generation.
     */
    public boolean isExtinct(Virus c) {
        int g1 = this.getTree(c).maxGen();
        int g2 = this.genOf(c);
        int g3 = this.maxGen();
        return g1 + g2 < g3;
    }

    /**
     * Recursively finds all Viruses which burst after time t and adds them
     * to a List.
     *
     * @param t the burstday after which we care about. Must be >0.
     * @return a List containing all Viruses born after time t.
     */
    public List<Virus> getAge(int t) {
        if (t < 0){
            throw new IllegalArgumentException("Burstday must be >=0.");
        }
        List<Virus> list = new LinkedList<>();
        if (this.founder.getBurstday() > t){
            list.add(this.founder);
        }
        for (VirusTree dt: children){
            list.addAll(dt.getAge(t));
        }
        return list;
    }

    /**
     * Sorts a list of VirusTrees by their burstdays, with the oldest VirusTree at index 0.
     *
     * @param list1 the List of VirusTrees to be sorted.
     * @return a List of VirusTrees sorted by age.
     */
    public static List<VirusTree> sortAge(List<VirusTree> list1){
        List<VirusTree> list2 = list1;
        Collections.sort(list2, new Comparator<VirusTree>() {
            @Override
            public int compare(VirusTree o1, VirusTree o2) {
                return Double.compare(o1.getFounder().getBurstday(), o2.getFounder().getBurstday());
            }
        });
        return list2;
    }
}
