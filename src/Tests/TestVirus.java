package Tests;

import hivSim.*;

import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.TestInstance;
import static org.junit.jupiter.api.Assertions.*;

@TestInstance(TestInstance.Lifecycle.PER_CLASS)
class TestVirus {
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
    public void testData() {
        d1223.setData("1223 MODIFIED");
        assertEquals("1222", d1222.getData());
        assertEquals("1223 MODIFIED", d1223.getData());
    }

    @Test
    public void testBurstday() {
        d1223.setBurstday(6);
        assertEquals(6, d1223.getBurstday());
        d12.setBurstday(0);
        assertEquals(0.0, d12.getBurstday());
        try {
            d13.setBurstday(-1);
        } catch (IllegalArgumentException e) {
            assertEquals(e.toString(), "java.lang.IllegalArgumentException: Burstday must be >= 0.0.");
        }
        d13.setBurstday(23.7);
        assertEquals(23.7, d13.getBurstday());
    }

    @Test
    public void testToString() {
        Virus testVirus = new Virus(null, 0);
        assertEquals("", testVirus.toString());
        assertEquals("1111", d1111.toString());
    }
}