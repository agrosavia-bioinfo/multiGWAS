package jmultigwas;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;

public class Model {

    RunAppWorker runAppWorker;

    Controller controller;

    public Model(Controller controller) {
        this.controller = controller;
    }

    public void runApplication(String outputDir, String inputValues) {
        System.out.println(">>> Running application...");
        System.out.println(">>> ouputDir: " + outputDir);
        String configFilenamePath = writeToFile(outputDir, inputValues);
        System.out.println(">>> config file: " + configFilenamePath);

        runAppWorker = new RunAppWorker(configFilenamePath, outputDir, controller);
        runAppWorker.execute();
    }

    String writeToFile(String outputDir, String inputValues) {
        String sep = File.separator;
        String dirName = new File(outputDir).getName();
        String configFilenamePath = outputDir + sep + dirName + ".config";

        try {
            FileWriter writer = new FileWriter(configFilenamePath);
            writer.write(inputValues);
            writer.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
        return (configFilenamePath);
    }
}
