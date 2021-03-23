package jmultigwas;

import java.awt.BorderLayout;
import java.awt.Desktop;
import java.awt.Dimension;
import java.io.File;
import java.io.IOException;
import java.net.URI;
import java.net.URISyntaxException;
import java.nio.file.Paths;
import javax.swing.JFrame;
import javax.swing.JMenu;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JTabbedPane;

class Controller extends JFrame {

    // Attributes
    Model model;
    JTabbedPane viewTabs;
    ViewInputs tabInputs;
    ViewFiles tabFiles;
    ViewToolBar viewToolBar;

    ViewOutputs tabOutputs;
    ViewResults tabResults;
    //JPanel tabInputs;

    JMenu menu, submenu;
    JMenuItem itemNew, itemOpen;

    // Methods
    public Controller(String text) {
        super(text);
        model = new Model(this);
        setTitle("JMultiGWAS Tool for GWAS");
        
    }

    public void init() {
        this.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        this.setSize(900, 500);

        initViewTabs();

        //this.setContentPane(viewTabs);
        this.add(viewTabs, BorderLayout.CENTER);
        this.add(viewToolBar, BorderLayout.WEST);

        // this.setMenu();
        this.setVisible(true);
        tabInputs.setEnabledInputs(false);
        //setTestMode (true);
    }
    
    // Start with test mode options
    public void setTestMode (boolean testMode) {
        if (testMode) {
            tabInputs.setTestMode (testMode);
            viewToolBar.setTestMode (testMode);
        }
    }

    public void initViewTabs() {
        viewTabs = new JTabbedPane();
        Dimension size = this.getSize();

        tabInputs = new ViewInputs(this, model);
        tabOutputs = new ViewOutputs(size);
        tabResults = new ViewResults(size);
        tabFiles = new ViewFiles(this);
        tabOutputs.init();

        viewTabs.addTab("Inputs", tabInputs);
        viewTabs.addTab("Outputs", tabOutputs);
        viewTabs.addTab("Results", tabResults);
        viewTabs.addTab("Files", tabFiles);
        
        viewToolBar = new ViewToolBar(this);
        
        //onNewProject();
    }

    public String getToolsToRun() {
        String tools = viewToolBar.getToolsToRun();
        return (tools);
    }

    public String getGeneAction() {
        String geneAction = viewToolBar.getGeneAction();
        return (geneAction);
    }

    public void onDefaultButton() {
        tabInputs.setDefaults();
    }

    public void onRunApplication() {
        if (tabInputs.checkCompleteInfo()==false) {
            JOptionPane.showMessageDialog(this, "Incomplete information", "MultiGWAS warning",
                    JOptionPane.WARNING_MESSAGE);
        }else {
            String outputDir  = tabInputs.getOutputDir();
            String values     = tabInputs.getInputValues();

            model.runApplication(outputDir, values);
            viewTabs.setSelectedIndex(1);
            tabOutputs.clearOutputs();
        }
    }

    public void onCancelButton() {
        viewTabs.setSelectedComponent(tabInputs);
    }

    public void onEndOfExecution() {
        System.out.println ("Hello");
        String SEP = File.separator;
        String workingDir = tabInputs.getOutputDir();
        String dirName = Paths.get(workingDir).getFileName().toString();
        String outputDir = workingDir + SEP + "out-" + dirName;
        
        // Get trait dir
        File[] directories = new File(outputDir).listFiles(File::isDirectory);
        String traitDir = directories[0].getName();
        System.out.println (">>>" + traitDir);
        
        String reportDir    = outputDir + SEP + traitDir + SEP + "report";
        String htmlFilename = outputDir + SEP + traitDir + SEP + "multiGWAS-report.html";
        //browseFile(htmlFilename, reportDir);
        
        writeLine ("Report of results in: <a href='file://"+htmlFilename+"'>"+htmlFilename+"</a>","html");
        
        // Set tabs 
        tabResults.showResults(htmlFilename);
        // viewTabs.setSelectedIndex(2);
        tabFiles.changeDir(reportDir);
        //tabResults = new ViewResults("/tmp/multiGWAS-report.html");  
    }
    public void onGenotypeFormat (String genoFormat) {
        if (genoFormat.matches("KMatrix|FitPoly|Updog")) {
            System.out.println ("Received KMatrix");
            tabInputs.enableMapComponents(true);
        }else if (genoFormat.matches("GWASpoly|VCF")) {
                tabInputs.enableMapComponents(false);
        }
    }

    public void browseFile(String url) {
        String myOS = System.getProperty("os.name").toLowerCase();
        OUT("(Your operating system is: " + myOS + ")\n");

        try {
            if (myOS.contains("win") && Desktop.isDesktopSupported()) { // Probably Windows
                OUT(" -- Going with Desktop.browse ...");
                Desktop desktop = Desktop.getDesktop();
                desktop.browse(new URI(url));
            } else { // Definitely Non-windows
                Runtime runtime = Runtime.getRuntime();
                if (myOS.contains("mac")) { // Apples
                    OUT(" -- Going on Apple with 'open'...");
                    runtime.exec("open " + url);
                } else if (myOS.contains("nix") || myOS.contains("nux")) { // Linux flavours 
                    OUT(" -- Going on Linux with 'xdg-open'...");
                    runtime.exec("xdg-open " + url);
                } else {
                    OUT("I was unable/unwilling to launch a browser in your OS :( #SadFace");
                }
            }
            OUT("\nThings have finished.");
        } catch (IOException | URISyntaxException eek) {
            OUT("**Stuff wrongly: " + eek.getMessage());
        }

    }

    public void writeLine(String s, String type) {
        tabOutputs.writeLine(s, type);
    }

    private void OUT(String string) {
        System.out.println(string);
    }

}
