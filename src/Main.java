import gurobi.GRBException;

import java.io.IOException;

public class Main {
    public static void main(String[] args) throws IOException, GRBException {
        Instance instance = readData.readData("loctrans.dat");
        LR lr = new LR(instance);
        lr.main();
        lr.report();
    }
}