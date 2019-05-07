package hivSim;

import org.neuroph.core.*;
import org.neuroph.core.data.*;
import org.neuroph.nnet.*;

import java.io.FileWriter;
import java.io.IOException;
import java.util.List;

public class NN {
    /**
     * Appends data (in calculate()'s return format) to a pre-existing .CSV
     * file.
     *
     * @param data an individual line of CSV data to be added.
     * @param fileName the name of the output file.
     */
    public static void appendData(List<Double> data, String fileName) throws IOException {
        FileWriter fw = new FileWriter(fileName + ".csv", true);
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < data.size() - 1; i++){
                sb.append(data.get(i).toString() + ",");
            }
            sb.append(data.get(data.size() - 1).toString() + '\n');
        fw.write(sb.toString());
        fw.close();
    }

    public static DataSet generateTrainingData(int M){
        return DataSet.createFromFile("C:\\Users\\Vruj\\Desktop\\hivSim\\data\\NN\\Train\\" + M + ".csv", M/2, 1,",", false);
    }

    public static DataSet generateTestingData(int M){
        return DataSet.createFromFile("C:\\Users\\Vruj\\Desktop\\hivSim\\data\\NN\\Test\\" + M + ".csv", M/2, 1, ",", false);
    }

    public static NeuralNetwork netTest(int M){
        NeuralNetwork net64 = new Perceptron(M/2, 1);
        DataSet trainingSet = generateTrainingData(M);
        net64.save("C:\\Users\\Vruj\\Desktop\\hivSim\\data\\NN\\net64.nnet");
        net64.learn(trainingSet);
        return net64;
    }


}
