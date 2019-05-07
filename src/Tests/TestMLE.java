package Tests;

import hivSim.MLE;

import java.io.FileNotFoundException;
import java.util.List;

import static org.junit.jupiter.api.Assertions.assertEquals;

public class TestMLE {

    public void testS_kmn(){
        MLE.S_kmn(1.01, 2.0, 5.0, 2.0);
    }

    public void testBinomial(){
        assertEquals(MLE.binomial(0,0), 1);
        assertEquals(MLE.binomial(4,2), 6);
        assertEquals(MLE.binomial(8,4), 70);
        assertEquals(MLE.binomial(0,0), 1);
    }
}
