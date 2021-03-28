/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package jmultigwas;

import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;
import java.util.concurrent.TimeUnit;
import javax.swing.SwingWorker;

/**
 *
 * @author lg
 */
public class RunAppWorker extends SwingWorker<Void, String> {

    Controller controller;
    String configFilenamePath;
    String outputDir;

    public RunAppWorker(String configFilenamePath, String outputDir, Controller controller) {
        this.controller = controller;
        this.configFilenamePath = configFilenamePath;
        this.outputDir = outputDir;
        
    }

    @Override
    protected Void doInBackground() throws Exception {
        ProcessBuilder processBuilder = new ProcessBuilder();
        processBuilder.directory(new File(outputDir));
        
        // -- Linux ... Run a shell command

        //processBuilder.command("bash", "-c", "pwd");
        //System.out.println ("....Out");
        //System.exit (0);
        Path outputPath       = Paths.get(outputDir);
        String outputDirName  = outputPath.getFileName().toString();
        String configFilename = outputDirName + ".config";
        
        String commandString = "multigwas.R " + configFilename;
        
        processBuilder.command("bash", "-c", commandString, outputPath.toString());
        System.out.println (">>> Command: " + processBuilder.command());
        
        try {
            Process p = processBuilder.start();

            try {
                BufferedInputStream bis = new BufferedInputStream(p.getInputStream());
                BufferedReader r = new BufferedReader(
                        new InputStreamReader(bis));

                String line;
                while ((line = r.readLine()) != null) {
                    System.out.println(line);
                    controller.writeLine (line,"");
                    if (line.contains("Moving")) {
                        //keywordFound = true;
                        break;
                    }
                }
            } finally {
                p.getInputStream().close();
            } 
        } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
        System.out.println ("END");
        controller.onEndOfExecution();
        
        return null;
    }
    

    @Override
    protected void process(List<String> chunks) {
        
    }

}
