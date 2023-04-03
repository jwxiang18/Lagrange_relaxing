import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;

public class readData {
    public static Instance readData(String filename) throws IOException {
        Instance instance = new Instance();
        try(BufferedReader br = Files.newBufferedReader(Paths.get(filename))){
            instance.buildlimit = Integer.parseInt(br.readLine());
            instance.ncities = Integer.parseInt(br.readLine());
            String[] line = br.readLine().trim().split("\\s+");
            instance.supply = new double[instance.ncities];
            for(int i = 0; i < instance.ncities ; i++){
                instance.supply[i] = Double.parseDouble(line[i]);
            }
            line = br.readLine().trim().split("\\s+");
            instance.demand = new double[instance.ncities];
            for (int i = 0; i < instance.ncities; i++) {
                instance.demand[i] = Double.parseDouble(line[i]);
            }
            instance.shipCost = new double[instance.ncities][instance.ncities];
            for (int i = 0; i < instance.ncities; i++) {
                line = br.readLine().trim().split("\\s+");
                for (int j = 0; j < instance.ncities; j++) {
                    instance.shipCost[i][j] = Double.parseDouble(line[j]);
                }
            }
            return instance;
        }
    }
}
