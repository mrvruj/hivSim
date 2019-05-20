package Tests;

import hivSim.*;
import nr.ran.*;

import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.TestInstance;

import static java.lang.StrictMath.log10;
import static org.junit.jupiter.api.Assertions.*;

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
//TODO: Finish test cases
@TestInstance(TestInstance.Lifecycle.PER_CLASS) //TODO What exactly does this do?
class TestVirusTree {
    private Virus<String> p1;
    private VirusTree Tree;

    private Virus<String> d11;
    private Virus<String> d12;
    private Virus<String> d13;

    private Virus<String> d111;
    private Virus<String> d112;
    private Virus<String> d113;

    private Virus<String> d121;
    private Virus<String> d122;

    private Virus<String> d1111;

    private Virus<String> d1131;
    private Virus<String> d1132;

    private Virus<String> d1211;
    private Virus<String> d1221;
    private Virus<String> d1222;
    private Virus<String> d1223;

    Ran rand = new Ran(314159267);

    /**
     * Generate a small viral phylogeny to test simulation methods.
     */
    @BeforeAll
    void setUp() {
        p1 = new Virus<>("1", 0);
        Tree = new VirusTree(p1);

        d11 = new Virus<>("11", 1);
        d12 = new Virus<>("12", 1);
        d13 = new Virus<>("13", 1);
        Tree.add(p1, d11);
        Tree.add(p1, d12);
        Tree.add(p1, d13);

        d111 = new Virus<>("111", 2);
        d112 = new Virus<>("112", 2);
        d113 = new Virus<>("113", 2);
        d121 = new Virus<>("121", 2);
        d122 = new Virus<>("122", 2);
        Tree.add(d11, d111);
        Tree.add(d11, d112);
        Tree.add(d11, d113);
        Tree.add(d12, d121);
        Tree.add(d12, d122);

        d1111 = new Virus<>("1111", 3);
        d1131 = new Virus<>("1131", 3);
        d1132 = new Virus<>("1132", 3);
        d1211 = new Virus<>("1211", 3);
        d1221 = new Virus<>("1221", 3);
        d1222 = new Virus<>("1222", 3);
        d1223 = new Virus<>("1223", 3);
        Tree.add(d111, d1111);
        Tree.add(d113, d1131);
        Tree.add(d113, d1132);
        Tree.add(d121, d1211);
        Tree.add(d122, d1221);
        Tree.add(d122, d1222);
        Tree.add(d122, d1223);
    }

    @Test
    public void testConstructors(){
        VirusTree TreeCopy = new VirusTree(Tree);
        assertEquals(3, Tree.getChildren().size());
        assertEquals(0, TreeCopy.getChildren().size()); //copies of a VirusTree share no children with their originals
        assertEquals(TreeCopy.getFounder(), Tree.getFounder());

        VirusTree TestTree = new VirusTree(p1); //a Virus can found multiple VirusTrees with mutually exclusive children
        assertEquals(0, TestTree.getChildren().size());
        assertNotEquals(0, Tree.getChildren().size());
    }

    @Test
    public void testGetFounder(){
        assertEquals(p1, Tree.getFounder());
        assertEquals(d13, Tree.getTree(d13).getFounder());
    }

    @Test
    public void testGetChildren(){
        assertTrue(Tree.getChildren().contains(Tree.getTree(d11)));
        assertTrue(Tree.getChildren().contains(Tree.getTree(d11)));
        assertTrue(Tree.getChildren().contains(Tree.getTree(d11)));
        assertEquals(0, Tree.getTree(d1111).getChildren().size());
    }

    @Test
    public void testGetTree(){
        assertEquals(d11, Tree.getTree(d11).getFounder());
        assertEquals(3, Tree.getTree(d11).getChildren().size());
        assertEquals(Tree, Tree.getTree(p1));
        assertNotEquals(Tree.getTree(d12), new VirusTree(d12));
    }

    @Test
    public void testContains() {
        assertTrue(Tree.getTree(d12).contains(d1221));
        assertFalse(Tree.getTree(d11).contains(d121));
        assertTrue(Tree.getTree(d111).contains(d1111));
        assertTrue(Tree.contains(p1));
        assertTrue(Tree.contains(d1111));
        Virus TestVirus = new Virus("test", 0.0);
        assertFalse(Tree.contains(TestVirus));
    }

    @Test
    public void testAdd(){
        Virus TestChild = new Virus("child", 314);
        Virus TestFounder = new Virus("founder", 302);
        VirusTree TestTree = new VirusTree(TestFounder);
        assertFalse(TestTree.contains(TestChild));
        TestTree.add(TestFounder, TestChild);
        assertTrue(TestTree.contains(TestChild));
        try{
            Tree.add(p1, null);
            Tree.add(null, TestChild);
        }
        catch(IllegalArgumentException e){
            assertEquals(e.toString(), "java.lang.IllegalArgumentException: Parent and child Viruses must be non-null.");
        }
    }

    @Test
    public void testSize(){
        assertEquals(16, Tree.size());
        assertEquals(7, Tree.getTree(d11).size());
        assertEquals(1, Tree.getTree(d1111).size());
    }

    @Test
    public void testToString(){
        assertEquals("1", Tree.toString());
        Tree.getFounder().setData("updated string");
        assertEquals("updated string", Tree.toString());
        Tree.getFounder().setData(null);
        assertEquals("", Tree.toString());
        Tree.getFounder().setData(1);   //reset for other tests
        assertEquals("1", Tree.toString());
    }

    @Test
    public void testRemoveChildren(){
        assertEquals(3, Tree.getTree(d122).getChildren().size());
        Tree.removeChildren(Tree.getTree(d122));
        assertEquals(0, Tree.getTree(d122).getChildren().size());
        Tree.removeChildren(Tree.getTree(d1111)); //remove all children from empty children List
        Tree.add(d122, d1221);
        Tree.add(d122, d1222);
        Tree.add(d122, d1223);
    }

    @Test
    public void testGetChildless() {
        List<VirusTree> ChildlessHolder = new ArrayList<>();
        ChildlessHolder.add(Tree.getTree(d1111));
        ChildlessHolder.add(Tree.getTree(d112));
        ChildlessHolder.add(Tree.getTree(d1131));
        ChildlessHolder.add(Tree.getTree(d1132));
        ChildlessHolder.add(Tree.getTree(d1211));
        ChildlessHolder.add(Tree.getTree(d1221));
        ChildlessHolder.add(Tree.getTree(d1222));
        ChildlessHolder.add(Tree.getTree(d1223));
        ChildlessHolder.add(Tree.getTree(d13));
        assertEquals(Tree.getChildless().size(), ChildlessHolder.size());
        for (VirusTree dt : Tree.getChildless()){
            assertTrue(ChildlessHolder.contains(dt));
        }
    }

    @Test
    public void testRouteTo() {
        LinkedList<Virus> test = new LinkedList<>();
        test.addFirst(d1111);
        test.addFirst(d111);
        test.addFirst(d11);
        test.addFirst(p1);
        assertEquals(Tree.routeTo(d1111), test);
        test.clear();
        test.add(p1);
        assertEquals(Tree.routeTo(p1), test);
    }

    @Test
    //TODO: test this more thoroughly
    public void testAddRandom(){
        Ran RNG = new Ran(1618033988);
        Poissondev poisson = new Poissondev(1.2, RNG.int64());
        Tree.addRandom(d1111, 4, true, 69.4, 34.7, poisson, RNG.int64());
        assertEquals(4, Tree.getTree(d1111).getChildren().size());
        VirusTree.removeChildren(Tree.getTree(d1111)); //this is unreachable if the line above fails, so testChildless() will likely also fail
    }

    @Test
    public void testSortAge(){
        List<VirusTree> sortThis = new ArrayList<>();
        sortThis.add(Tree.getTree(d1111));
        sortThis.add(Tree.getTree(d11));
        sortThis.add(Tree.getTree(d121));
        sortThis.add(Tree.getTree(d1223));
        sortThis.add(Tree.getTree(d13));
        sortThis.add(Tree.getTree(p1));
        List<VirusTree> sorted = new ArrayList<>();
        sorted.add(Tree.getTree(p1));
        sorted.add(Tree.getTree(d11));
        sorted.add(Tree.getTree(d13));
        sorted.add(Tree.getTree(d121));
        sorted.add(Tree.getTree(d1111));
        sorted.add(Tree.getTree(d1223));
        assertEquals(sorted.toString(), VirusTree.sortAge(sortThis).toString());
    }

    @Test
    public void testPhylogenySim(){
        double r0 = 2.5937424601; double mu = 0.0551; int M = 1024; int maxVirus = 6000; boolean cTildas = true;
        long randomLong = 271828182845904523L;
        //long randomLong = System.currentTimeMillis();
        Ran RNG = new Ran(randomLong);
        CSV.setTolerance(6);
        MLE.populateBinom(1024);
        MLE.populateSkmn(2, M, M);

        Poissondev poisson = new Poissondev(r0, randomLong);
        VirusTree FoundingVirusTree = new VirusTree(new Virus('0', 0));

        List<VirusTree> live = VirusTree.phylogenySim(FoundingVirusTree, maxVirus, r0, true, randomLong, poisson);
        List<Virus> sample = VirusTree.sample(live, M, RNG.int64());
        List<Integer> eta = VirusTree.siteFrequencySpectrum(VirusTree.ancestorFrequencySpectrum(sample, FoundingVirusTree), mu, poisson);
        List<Integer> etaTilda = VirusTree.foldedSiteFrequencySpectrum(eta);

        double r_MLE = MLE.rHat(etaTilda, M, mu, cTildas); //can we get rid of M dependency?
        double s_hat_MLE = MLE.sigmaHat(r0, M, mu, etaTilda, cTildas);

        double r_MLE_EM = MLE.rHatEM(M, mu, etaTilda);
        double s_hat_MLE_EM = MLE.sigmaHatEM(r_MLE_EM, M, etaTilda);

        double r_MOM = MOM.rHat(mu, M,etaTilda);
        double s_hat_MOM = MOM.sigmaHat(M, r_MOM,mu);

        double log10e = log10(2.718281828459045235360287471);
        System.out.printf("%n etaTilda: %n%s%n", etaTilda.toString());

        System.out.printf("%nM: %s%n", M);
        System.out.printf("r0: %f%n", r0);
        System.out.printf("%nMLE r_hat: %s%n", r_MLE);
        System.out.printf("MLE s_hat: %s%n%n", log10e*s_hat_MLE);

        System.out.printf("MLE_EM r_hat: %s%n", r_MLE_EM);
        System.out.printf("MLE_EM s_hat: %s%n%n", log10e*s_hat_MLE_EM);

        System.out.printf("MOM r_hat: %s%n", r_MOM);
        System.out.printf("MOM s_hat: %s%n%n", log10e*s_hat_MOM);
    }

    @Test
    public void testBinomial(){
        CSV.setTolerance(9);
        MLE.populateBinom(512);
        System.out.println(MLE.binomial(512, 256));
    }

    @Test
    public void testATilda(){
        CSV.setTolerance(9);
        MLE.populateBinom(1024);
        MLE.populateSkmn(2, 1024, 1024);
        double x = MLE.aTilda(1.21, 512, 2);
        System.out.println(x);
    }

    @Test
    public void testIsOlderThan(){
        assertTrue(Tree.getTree(d111).isOlderThan(Tree.getTree(d1111)));
        assertFalse(Tree.getTree(d111).isOlderThan(Tree.getTree(d11)));
        assertTrue(Tree.isOlderThan(Tree.getTree(d111)));
        assertTrue(Tree.isOlderThan(Tree));
        assertTrue(Tree.getTree(d11).isOlderThan(Tree.getTree(d12)));
    }

    @Test
    public void testFactorial(){
        long old = System.currentTimeMillis();
        double g = MLE.factorialG(1024);
        System.out.printf("Gamma time: %f%n", (System.currentTimeMillis() - old)/1000.0);
        long new1 = System.currentTimeMillis();
        double b = MLE.factorial(1024).doubleValue();
        System.out.printf("BigI time: %f%n", (System.currentTimeMillis() - new1)/1000.0);
        double d = g - b;
        System.out.printf("Gamma: %f; BigInteger: %f; diff = %f", g, b, d);
    }
    //TODO: finish testing these


    //@Test
    public void testSample() {
    }

    //@Test
    public void testAncestryCounter(){
    }

    //@Test
    public void testAncestorFrequencySpectrum(){
    }

    //@Test
    public void testSampleFrequencySpectrum(){
    }

    //@Test
    public void testFoldedSampleFrequencySpectrum(){
    }

    //Unused functions below

    public void testGenOf() {
        assertEquals(Tree.genOf(p1), 0);
        assertEquals(Tree.genOf(d112), 2);
        assertEquals(Tree.genOf(d1221), 3);
    }

    public void testMaxGen() {
        assertEquals(Tree.maxGen(), 3);
        assertEquals(Tree.getTree(d11).maxGen(), 2);
    }

    public void testWidthAtDepth(){
    }

    public void testIsExtinct() {
        assertEquals(Tree.isExtinct(d112), true);
        assertEquals(Tree.isExtinct(d11), false);
        assertEquals(Tree.isExtinct(p1), false);
        assertEquals(Tree.isExtinct(d113), false);
    }

    public void testSharedAncestor() {
        assertEquals(Tree.sharedAncestorOf(d11, d12), p1);
        assertEquals(Tree.sharedAncestorOf(d1132, d1222), p1);
        assertEquals(Tree.sharedAncestorOf(d111, d113), d11);
    }

    public void testGetGen() {
        System.out.println(Tree.getGen(0).toString());
        System.out.println(Tree.getGen(1).toString());
        System.out.println(Tree.getGen(2).toString());
        System.out.println(Tree.getGen(Tree.maxGen()).toString());
    }

    public void testGetAge() {
        System.out.println(Tree.getAge(0).toString());
        System.out.println(Tree.getAge(1).toString());
        System.out.println(Tree.getAge(2).toString());
        System.out.println(Tree.getAge(Tree.maxGen()).toString());
    }

    public void testGetAll(){
        System.out.println(Tree.getAll(Tree).toString());
    }

    public void testGetParent(){
        assertEquals(p1, Tree.getParent(d11));
        assertEquals(Tree.getParent(d1211), d121);
        assertEquals(Tree.getParent(p1), null);
    }
}
