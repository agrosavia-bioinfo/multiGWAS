package jmultigwas;

import java.awt.Component;
import java.io.File;
import java.util.prefs.Preferences;
import javax.swing.JFileChooser;
import javax.swing.JOptionPane;

public class ViewInputs extends javax.swing.JPanel {

    // Attributes
    Model model;
    Controller controller;
    String outputDir;
    Preferences prefs;

    final JFileChooser fc = new JFileChooser();
    String LAST_USED_FOLDER = "";

    // Methods
    public ViewInputs(Controller controller, Model model) {
        this.controller = controller;
        this.model = model;
        initComponents();
        setDefaults();
        fieldOutputDir.setText(getLastDir().toString());
        //fieldGeno.setText(getLastDir().toString());
        //fieldPheno.setText(getLastDir().toString());

    }

    public File getLastDir() {
        prefs = Preferences.userRoot().node(getClass().getName());
        File curFile = new File(prefs.get(LAST_USED_FOLDER, new File(".").getAbsolutePath()));
        return curFile;
    }

    public void setDefaults() {
        fieldOutputDir.setText(System.getProperty("user.home"));
        fieldGeno.setText("");
        fieldPheno.setText("");
        fieldMap.setText("");

        fieldPloidy.setSelectedIndex(0);
        fieldSignificance.setText("0.05");
        fieldCorrection.setSelectedIndex(0);
        fieldModel.setSelectedIndex(0);
        fieldSNPs.setSelectedIndex(9);

        fieldFiltering.setSelectedIndex(0);
        fieldFilterMAF.setText("0.01");
        fieldFilterMIND.setText("0.1");
        fieldFilterGENO.setText("0.1");
        fieldFilterHWE.setText("1e-10");
        fieldR2.setText ("0.9");
    }

    public void clearInputs() {
        fieldOutputDir.setText("");
        fieldGeno.setText("");
        fieldPheno.setText("");
        fieldMap.setText("");

        fieldPloidy.setSelectedIndex(0);
        fieldSignificance.setText("");
        fieldCorrection.setSelectedIndex(0);
        fieldModel.setSelectedIndex(0);
        fieldSNPs.setSelectedIndex(9);

        fieldFiltering.setSelectedIndex(0);
        fieldFilterMAF.setText("");
        fieldFilterMIND.setText("");
        fieldFilterGENO.setText("");
        fieldFilterHWE.setText("");
        fieldR2.setText("");
        setEnabledInputs(false);
    }

    public void setEnabledInputs(Boolean flag) {
        Component components[] = this.getComponents();

        for (Component c : panelPaths.getComponents()) {
            c.setEnabled(flag);
        }
        enableMapComponents (false);

        for (Component c : panelParameters.getComponents()) {
            c.setEnabled(flag);
        }

        for (Component c : panelFilters.getComponents()) {
            c.setEnabled(flag);
        }

        fieldOutputDir.setEnabled(true);
        labelOutputDir.setEnabled(true);
        selOutputDirButton.setEnabled(true);     
    }
    public void enableMapComponents (boolean enableFlag) {
        mapLabel.setEnabled (enableFlag);
        fieldMap.setEnabled (enableFlag);
        mapSelButton.setEnabled (enableFlag);
    } 
    
    // Init form with test mode options
    public void setTestMode (boolean testMode) {
        if (testMode) {
            fieldGeno.setText ("/home/lg/tpm/test/gwaspoly-genotype-500.csv");
            fieldPheno.setText ("/home/lg/tmp/test/example-phenotype-single-trait.csv");
            fieldFiltering.setSelectedIndex(1);
            
        }
    }

    public String getInputValues() {
        StringBuffer txt = new StringBuffer(500);
        String ln = System.lineSeparator();
        txt.append("# Files:" + ln);
        txt.append("  genotypeFile        : " + '"' + fieldGeno.getText() + '"' + ln);
        txt.append("  phenotypeFile       : " + '"' + fieldPheno.getText() + '"' + ln);
        txt.append("  genotypeFormat      : " + '"' + genotypeFormatCBox.getSelectedItem().toString() + '"' + ln);        
        txt.append("  mapFile             : " + '"' + fieldMap.getText() + '"' + ln); 

        txt.append("# GWAS model:" + ln);
        txt.append("  ploidy              : " + fieldPloidy.getSelectedItem().toString() + ln);        
        txt.append("  significanceLevel   : " + fieldSignificance.getText() + ln);
        txt.append("  correctionMethod    : " + '"' + fieldCorrection.getSelectedItem().toString() + '"' + ln);
        txt.append("  gwasModel           : " + fieldModel.getSelectedItem().toString() + ln);
        txt.append("  geneAction          : " + '"' + controller.getGeneAction() + '"' + ln);
        txt.append("  R2                  : " + '"' + fieldR2.getText() + '"' + ln);
        
        txt.append("# Quality control:" + ln);
        txt.append("  filtering           : " + fieldFiltering.getSelectedItem().toString() + ln);
        txt.append("  MAF                 : " + fieldFilterMAF.getText() + ln);
        txt.append("  MIND                : " + fieldFilterMIND.getText() + ln);
        txt.append("  GENO                : " + fieldFilterGENO.getText() + ln);
        txt.append("  HWE                 : " + fieldFilterHWE.getText() + ln);

        txt.append("# Tools:" + ln);
        txt.append("  nBest               : " + fieldSNPs.getSelectedItem().toString() + ln);
        txt.append("  tools               : " + '"' + controller.getToolsToRun() + '"' + ln);

        return txt.toString();
    }

    public boolean checkCompleteInfo() {
        if (fieldOutputDir.getText().equals ("")) return false;
        if (fieldGeno.getText().equals ("")) return false;
        if (fieldPheno.getText().equals ("")) return false;
        if (fieldSignificance.getText().equals ("")) return false;
        if (fieldFilterMAF.getText().equals ("")) return false;
        if (fieldFilterMIND.getText().equals ("")) return false;
        if (fieldFilterGENO.getText().equals ("")) return false;
        if (fieldFilterHWE.getText().equals ("")) return false;
        
        if (genotypeFormatCBox.getSelectedItem().toString().equals ("KMatrix") ||
            genotypeFormatCBox.getSelectedItem().toString().equals ("FitPoly") ||
            genotypeFormatCBox.getSelectedItem().toString().equals ("Updog"))
                if (fieldMap.getText().equals("")) return false;
                
        if (controller.getToolsToRun().equals ("")) return false;
        System.out.println("Success: complete information!!");
        return true;
    }

    public void runApplication() {
        String outputDir = fieldOutputDir.getText();
        String fieldsText = getInputValues();

        //controller.onRunApplication(outputDir, fieldsText);

        //model.runApplication(outputDir, fieldsText.toString());
    }

    public String getOutputDir() {
        outputDir = fieldOutputDir.getText();
        return outputDir;
    }

    /**
     * This method is called from within the constructor to initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is always
     * regenerated by the Form Editor.
     */
    @SuppressWarnings("unchecked")
    // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
    private void initComponents() {

        panelInputs = new javax.swing.JPanel();
        panelFilesTitle = new javax.swing.JPanel();
        labelInputOutput = new javax.swing.JLabel();
        panelPaths = new javax.swing.JPanel();
        labelOutputDir = new javax.swing.JLabel();
        fieldOutputDir = new javax.swing.JTextField();
        fieldPheno = new javax.swing.JTextField();
        jLabel2 = new javax.swing.JLabel();
        genotypeLabel = new javax.swing.JLabel();
        fieldGeno = new javax.swing.JTextField();
        selGenotypeButton = new javax.swing.JButton();
        selOutputDirButton = new javax.swing.JButton();
        selPhenotypeButton = new javax.swing.JButton();
        panelGenotypeInfo = new javax.swing.JInternalFrame();
        jCheckBox1 = new javax.swing.JCheckBox();
        jLabel9 = new javax.swing.JLabel();
        jTextField1 = new javax.swing.JTextField();
        genotypeFormatLabel = new javax.swing.JLabel();
        genotypeFormatCBox = new javax.swing.JComboBox<>();
        mapLabel = new javax.swing.JLabel();
        fieldMap = new javax.swing.JTextField();
        mapSelButton = new javax.swing.JButton();
        panelOptions = new javax.swing.JPanel();
        panelParametersTitle = new javax.swing.JPanel();
        jLabel1 = new javax.swing.JLabel();
        panelParameters = new javax.swing.JPanel();
        labelGwasModel = new javax.swing.JLabel();
        fieldModel = new javax.swing.JComboBox<>();
        jLabel3 = new javax.swing.JLabel();
        fieldSignificance = new javax.swing.JTextField();
        jLabel4 = new javax.swing.JLabel();
        fieldCorrection = new javax.swing.JComboBox<>();
        jLabel5 = new javax.swing.JLabel();
        fieldSNPs = new javax.swing.JComboBox<>();
        jLabel6 = new javax.swing.JLabel();
        fieldFiltering = new javax.swing.JComboBox<>();
        labelPloidy = new javax.swing.JLabel();
        fieldPloidy = new javax.swing.JComboBox<>();
        labelR2 = new javax.swing.JLabel();
        fieldR2 = new javax.swing.JTextField();
        filler1 = new javax.swing.Box.Filler(new java.awt.Dimension(0, 0), new java.awt.Dimension(15, 0), new java.awt.Dimension(32767, 0));
        panelFiltersTitle = new javax.swing.JPanel();
        panelFilters = new javax.swing.JPanel();
        fieldFilterMIND = new javax.swing.JTextField();
        fieldFilterGENO = new javax.swing.JTextField();
        fieldFilterHWE = new javax.swing.JTextField();
        fieldFilterMAF = new javax.swing.JTextField();
        jLabel7 = new javax.swing.JLabel();
        jLabel15 = new javax.swing.JLabel();
        jLabel16 = new javax.swing.JLabel();
        jLabel17 = new javax.swing.JLabel();
        jLabel18 = new javax.swing.JLabel();

        setBackground(new java.awt.Color(153, 153, 255));
        setPreferredSize(new java.awt.Dimension(780, 650));
        setLayout(new java.awt.BorderLayout());

        panelInputs.setPreferredSize(new java.awt.Dimension(780, 650));
        panelInputs.setLayout(new java.awt.BorderLayout());

        panelFilesTitle.setBorder(javax.swing.BorderFactory.createLineBorder(new java.awt.Color(0, 0, 0)));
        panelFilesTitle.setLayout(new java.awt.BorderLayout());

        labelInputOutput.setBackground(javax.swing.UIManager.getDefaults().getColor("Button.disabledToolBarBorderBackground"));
        labelInputOutput.setText("Input/Output:");
        labelInputOutput.setOpaque(true);
        panelFilesTitle.add(labelInputOutput, java.awt.BorderLayout.PAGE_START);

        panelPaths.setBorder(javax.swing.BorderFactory.createLineBorder(new java.awt.Color(0, 0, 0)));

        labelOutputDir.setText("Working dir:");

        fieldOutputDir.setText("/home/lg/AAA");
        fieldOutputDir.setToolTipText("Select genotype file");

        fieldPheno.setToolTipText("Select phenotype file");
        fieldPheno.setPreferredSize(new java.awt.Dimension(90, 19));

        jLabel2.setText("Phenotype file:");

        genotypeLabel.setText("Genotype file:");

        fieldGeno.setToolTipText("Select genotype file");
        fieldGeno.setPreferredSize(new java.awt.Dimension(90, 19));

        selGenotypeButton.setText("Select...");
        selGenotypeButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                selGenotypeButtonActionPerformed(evt);
            }
        });

        selOutputDirButton.setText("Select...");
        selOutputDirButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                selOutputDirButtonActionPerformed(evt);
            }
        });

        selPhenotypeButton.setText("Select...");
        selPhenotypeButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                selPhenotypeButtonActionPerformed(evt);
            }
        });

        panelGenotypeInfo.setBorder(javax.swing.BorderFactory.createLineBorder(new java.awt.Color(0, 0, 0)));
        panelGenotypeInfo.setTitle("Genotype information:");
        panelGenotypeInfo.setVisible(true);

        jCheckBox1.setText("Non-model organism");

        jLabel9.setText("Number of chromosomes:");

        jTextField1.setText("1");

        genotypeFormatLabel.setText("Genotype format:");

        genotypeFormatCBox.setModel(new javax.swing.DefaultComboBoxModel<>(new String[] { "GWASpoly", "KMatrix", "VCF", "FitPoly", "Updog", " " }));
        genotypeFormatCBox.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                genotypeFormatCBoxActionPerformed(evt);
            }
        });

        mapLabel.setText("Map file:");
        mapLabel.setEnabled(false);

        fieldMap.setToolTipText("Select phenotype file");
        fieldMap.setEnabled(false);
        fieldMap.setPreferredSize(new java.awt.Dimension(90, 19));

        mapSelButton.setText("Select...");
        mapSelButton.setEnabled(false);
        mapSelButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                mapSelButtonActionPerformed(evt);
            }
        });

        javax.swing.GroupLayout panelGenotypeInfoLayout = new javax.swing.GroupLayout(panelGenotypeInfo.getContentPane());
        panelGenotypeInfo.getContentPane().setLayout(panelGenotypeInfoLayout);
        panelGenotypeInfoLayout.setHorizontalGroup(
            panelGenotypeInfoLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(panelGenotypeInfoLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(panelGenotypeInfoLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(panelGenotypeInfoLayout.createSequentialGroup()
                        .addComponent(genotypeFormatLabel)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                        .addComponent(genotypeFormatCBox, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addGap(32, 32, 32)
                        .addComponent(jCheckBox1)
                        .addGap(18, 18, 18)
                        .addComponent(jLabel9)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(jTextField1, javax.swing.GroupLayout.PREFERRED_SIZE, 36, javax.swing.GroupLayout.PREFERRED_SIZE))
                    .addGroup(panelGenotypeInfoLayout.createSequentialGroup()
                        .addComponent(mapLabel)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                        .addComponent(fieldMap, javax.swing.GroupLayout.PREFERRED_SIZE, 447, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(mapSelButton)))
                .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );
        panelGenotypeInfoLayout.setVerticalGroup(
            panelGenotypeInfoLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(panelGenotypeInfoLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(panelGenotypeInfoLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(genotypeFormatLabel)
                    .addComponent(genotypeFormatCBox, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(jCheckBox1)
                    .addComponent(jLabel9)
                    .addComponent(jTextField1, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addGroup(panelGenotypeInfoLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(mapLabel)
                    .addComponent(fieldMap, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(mapSelButton, javax.swing.GroupLayout.PREFERRED_SIZE, 19, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );

        javax.swing.GroupLayout panelPathsLayout = new javax.swing.GroupLayout(panelPaths);
        panelPaths.setLayout(panelPathsLayout);
        panelPathsLayout.setHorizontalGroup(
            panelPathsLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(panelPathsLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(panelPathsLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(panelPathsLayout.createSequentialGroup()
                        .addComponent(labelOutputDir)
                        .addGap(47, 47, 47)
                        .addComponent(fieldOutputDir)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(selOutputDirButton))
                    .addGroup(panelPathsLayout.createSequentialGroup()
                        .addGroup(panelPathsLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(genotypeLabel)
                            .addComponent(jLabel2))
                        .addGap(30, 30, 30)
                        .addGroup(panelPathsLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(fieldGeno, javax.swing.GroupLayout.DEFAULT_SIZE, 515, Short.MAX_VALUE)
                            .addComponent(fieldPheno, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addGroup(panelPathsLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(selPhenotypeButton, javax.swing.GroupLayout.Alignment.TRAILING)
                            .addComponent(selGenotypeButton, javax.swing.GroupLayout.Alignment.TRAILING)))
                    .addComponent(panelGenotypeInfo)))
        );
        panelPathsLayout.setVerticalGroup(
            panelPathsLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(panelPathsLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(panelPathsLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(selOutputDirButton, javax.swing.GroupLayout.PREFERRED_SIZE, 19, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(labelOutputDir)
                    .addComponent(fieldOutputDir, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addGap(12, 12, 12)
                .addGroup(panelPathsLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(genotypeLabel)
                    .addComponent(fieldGeno, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(selGenotypeButton, javax.swing.GroupLayout.PREFERRED_SIZE, 19, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addGroup(panelPathsLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(jLabel2)
                    .addComponent(fieldPheno, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(selPhenotypeButton, javax.swing.GroupLayout.PREFERRED_SIZE, 19, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addComponent(panelGenotypeInfo, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );

        panelFilesTitle.add(panelPaths, java.awt.BorderLayout.CENTER);

        panelInputs.add(panelFilesTitle, java.awt.BorderLayout.PAGE_START);

        panelOptions.setBorder(new javax.swing.border.MatteBorder(null));
        panelOptions.setLayout(null);

        panelParametersTitle.setBorder(javax.swing.BorderFactory.createLineBorder(new java.awt.Color(0, 0, 0)));
        panelParametersTitle.setPreferredSize(new java.awt.Dimension(352, 265));
        panelParametersTitle.setLayout(new java.awt.BorderLayout());

        jLabel1.setBackground(javax.swing.UIManager.getDefaults().getColor("Button.disabledToolBarBorderBackground"));
        jLabel1.setText("GWAS Parameters:");
        jLabel1.setOpaque(true);
        panelParametersTitle.add(jLabel1, java.awt.BorderLayout.PAGE_START);

        panelParameters.setBorder(javax.swing.BorderFactory.createLineBorder(new java.awt.Color(0, 0, 0)));
        panelParameters.setAlignmentY(1.5F);
        panelParameters.setPreferredSize(new java.awt.Dimension(350, 250));

        labelGwasModel.setText("GWAS model:");

        fieldModel.setModel(new javax.swing.DefaultComboBoxModel<>(new String[] { "Full", "Naive" }));

        jLabel3.setText("Significance Level:");

        fieldSignificance.setText("0.05");

        jLabel4.setText("Correction Method:");

        fieldCorrection.setModel(new javax.swing.DefaultComboBoxModel<>(new String[] { "Bonferroni", "FDR" }));
        fieldCorrection.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                fieldCorrectionActionPerformed(evt);
            }
        });

        jLabel5.setText("Number of Best SNPs:");

        fieldSNPs.setModel(new javax.swing.DefaultComboBoxModel<>(new String[] { "1", "2", "3", "4", "5", "6", "7", "8", "9", "10" }));
        fieldSNPs.setSelectedIndex(9);

        jLabel6.setText("Filtering:");

        fieldFiltering.setModel(new javax.swing.DefaultComboBoxModel<>(new String[] { "TRUE", "FALSE" }));
        fieldFiltering.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                fieldFilteringActionPerformed(evt);
            }
        });

        labelPloidy.setText("Ploidy:");

        fieldPloidy.setModel(new javax.swing.DefaultComboBoxModel<>(new String[] { "4", "2" }));

        labelR2.setText("LD threshold (R2):");

        fieldR2.setText("1.0");
        fieldR2.setPreferredSize(new java.awt.Dimension(30, 18));
        fieldR2.setRequestFocusEnabled(false);

        javax.swing.GroupLayout panelParametersLayout = new javax.swing.GroupLayout(panelParameters);
        panelParameters.setLayout(panelParametersLayout);
        panelParametersLayout.setHorizontalGroup(
            panelParametersLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(panelParametersLayout.createSequentialGroup()
                .addGap(9, 9, 9)
                .addGroup(panelParametersLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(panelParametersLayout.createSequentialGroup()
                        .addComponent(labelPloidy)
                        .addGap(122, 122, 122)
                        .addComponent(fieldPloidy, javax.swing.GroupLayout.PREFERRED_SIZE, 130, javax.swing.GroupLayout.PREFERRED_SIZE))
                    .addGroup(panelParametersLayout.createSequentialGroup()
                        .addComponent(labelGwasModel)
                        .addGap(76, 76, 76)
                        .addComponent(fieldModel, javax.swing.GroupLayout.PREFERRED_SIZE, 130, javax.swing.GroupLayout.PREFERRED_SIZE))
                    .addGroup(panelParametersLayout.createSequentialGroup()
                        .addComponent(jLabel3)
                        .addGap(39, 39, 39)
                        .addComponent(fieldSignificance, javax.swing.GroupLayout.PREFERRED_SIZE, 130, javax.swing.GroupLayout.PREFERRED_SIZE))
                    .addGroup(panelParametersLayout.createSequentialGroup()
                        .addComponent(jLabel4)
                        .addGap(33, 33, 33)
                        .addComponent(fieldCorrection, javax.swing.GroupLayout.PREFERRED_SIZE, 130, javax.swing.GroupLayout.PREFERRED_SIZE))
                    .addGroup(panelParametersLayout.createSequentialGroup()
                        .addComponent(jLabel5)
                        .addGap(16, 16, 16)
                        .addComponent(fieldSNPs, javax.swing.GroupLayout.PREFERRED_SIZE, 130, javax.swing.GroupLayout.PREFERRED_SIZE))
                    .addGroup(panelParametersLayout.createSequentialGroup()
                        .addComponent(jLabel6)
                        .addGap(106, 106, 106)
                        .addComponent(fieldFiltering, javax.swing.GroupLayout.PREFERRED_SIZE, 130, javax.swing.GroupLayout.PREFERRED_SIZE))
                    .addGroup(panelParametersLayout.createSequentialGroup()
                        .addComponent(labelR2, javax.swing.GroupLayout.PREFERRED_SIZE, 130, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addGap(40, 40, 40)
                        .addComponent(fieldR2, javax.swing.GroupLayout.PREFERRED_SIZE, 130, javax.swing.GroupLayout.PREFERRED_SIZE))))
        );
        panelParametersLayout.setVerticalGroup(
            panelParametersLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(panelParametersLayout.createSequentialGroup()
                .addGap(19, 19, 19)
                .addGroup(panelParametersLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(labelPloidy)
                    .addComponent(fieldPloidy, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addGap(6, 6, 6)
                .addGroup(panelParametersLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(labelGwasModel)
                    .addComponent(fieldModel, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addGap(6, 6, 6)
                .addGroup(panelParametersLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(jLabel3)
                    .addComponent(fieldSignificance, javax.swing.GroupLayout.PREFERRED_SIZE, 20, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addGap(10, 10, 10)
                .addGroup(panelParametersLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(jLabel4)
                    .addComponent(fieldCorrection, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addGap(6, 6, 6)
                .addGroup(panelParametersLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(jLabel5)
                    .addComponent(fieldSNPs, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addGap(6, 6, 6)
                .addGroup(panelParametersLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(jLabel6)
                    .addComponent(fieldFiltering, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addGap(6, 6, 6)
                .addGroup(panelParametersLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(labelR2)
                    .addComponent(fieldR2, javax.swing.GroupLayout.PREFERRED_SIZE, 20, javax.swing.GroupLayout.PREFERRED_SIZE)))
        );

        panelParametersTitle.add(panelParameters, java.awt.BorderLayout.CENTER);

        panelOptions.add(panelParametersTitle);
        panelParametersTitle.setBounds(13, 13, 352, 265);
        panelOptions.add(filler1);
        filler1.setBounds(405, 138, 15, 0);

        panelFiltersTitle.setBorder(javax.swing.BorderFactory.createLineBorder(new java.awt.Color(0, 0, 0)));
        panelFiltersTitle.setPreferredSize(new java.awt.Dimension(352, 265));
        panelFiltersTitle.setLayout(new java.awt.BorderLayout());

        panelFilters.setBorder(javax.swing.BorderFactory.createLineBorder(new java.awt.Color(0, 0, 0)));

        fieldFilterMIND.setText("0.1");
        fieldFilterMIND.setPreferredSize(new java.awt.Dimension(40, 15));

        fieldFilterGENO.setText("0.1");
        fieldFilterGENO.setPreferredSize(new java.awt.Dimension(40, 15));
        fieldFilterGENO.setRequestFocusEnabled(false);

        fieldFilterHWE.setText("1e-10");
        fieldFilterHWE.setPreferredSize(new java.awt.Dimension(40, 15));

        fieldFilterMAF.setText("0.01");
        fieldFilterMAF.setPreferredSize(new java.awt.Dimension(25, 15));
        fieldFilterMAF.setRequestFocusEnabled(false);
        fieldFilterMAF.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                fieldFilterMAFActionPerformed(evt);
            }
        });

        jLabel7.setText("<html>Minimum Minor Allele Frequency (MAF) for a SNP to be kept: </html>");
        jLabel7.setName("<html>Exclude SNPs with minor <br> allele frequency less than </html>"); // NOI18N

        jLabel15.setText("<html>Maximum proportion of missing values for a SNP to be kept </html>");
        jLabel15.setName("<html>Exclude SNPs with minor <br> allele frequency less than </html>"); // NOI18N

        jLabel16.setText("<html>Maximum proportion of missing values for a sample to be kept: </html>");
        jLabel16.setName("<html>Exclude SNPs with minor <br> allele frequency less than </html>"); // NOI18N
        jLabel16.setPreferredSize(new java.awt.Dimension(440, 15));

        jLabel17.setText("<html>Filters out SNPs with HWE exact test p-value below threshold: </html>");
        jLabel17.setName("<html>Exclude SNPs with minor <br> allele frequency less than </html>"); // NOI18N

        javax.swing.GroupLayout panelFiltersLayout = new javax.swing.GroupLayout(panelFilters);
        panelFilters.setLayout(panelFiltersLayout);
        panelFiltersLayout.setHorizontalGroup(
            panelFiltersLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(panelFiltersLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(panelFiltersLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(panelFiltersLayout.createSequentialGroup()
                        .addComponent(jLabel16, javax.swing.GroupLayout.PREFERRED_SIZE, 232, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(fieldFilterGENO, javax.swing.GroupLayout.PREFERRED_SIZE, 80, javax.swing.GroupLayout.PREFERRED_SIZE))
                    .addGroup(panelFiltersLayout.createSequentialGroup()
                        .addComponent(jLabel17, javax.swing.GroupLayout.PREFERRED_SIZE, 232, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(fieldFilterHWE, javax.swing.GroupLayout.PREFERRED_SIZE, 80, javax.swing.GroupLayout.PREFERRED_SIZE))
                    .addGroup(panelFiltersLayout.createSequentialGroup()
                        .addComponent(jLabel15, javax.swing.GroupLayout.PREFERRED_SIZE, 232, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(fieldFilterMIND, javax.swing.GroupLayout.PREFERRED_SIZE, 80, javax.swing.GroupLayout.PREFERRED_SIZE))
                    .addGroup(panelFiltersLayout.createSequentialGroup()
                        .addComponent(jLabel7, javax.swing.GroupLayout.PREFERRED_SIZE, 232, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(fieldFilterMAF, javax.swing.GroupLayout.PREFERRED_SIZE, 80, javax.swing.GroupLayout.PREFERRED_SIZE)))
                .addContainerGap(50, Short.MAX_VALUE))
        );
        panelFiltersLayout.setVerticalGroup(
            panelFiltersLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(panelFiltersLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(panelFiltersLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(jLabel7, javax.swing.GroupLayout.PREFERRED_SIZE, 40, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(fieldFilterMAF, javax.swing.GroupLayout.PREFERRED_SIZE, 40, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addGroup(panelFiltersLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                    .addComponent(jLabel15)
                    .addComponent(fieldFilterMIND, javax.swing.GroupLayout.PREFERRED_SIZE, 40, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addGroup(panelFiltersLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                    .addComponent(jLabel16, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addComponent(fieldFilterGENO, javax.swing.GroupLayout.PREFERRED_SIZE, 40, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addGroup(panelFiltersLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                    .addComponent(fieldFilterHWE, javax.swing.GroupLayout.PREFERRED_SIZE, 40, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(jLabel17))
                .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );

        jLabel7.getAccessibleContext().setAccessibleName("");

        panelFiltersTitle.add(panelFilters, java.awt.BorderLayout.CENTER);

        jLabel18.setBackground(javax.swing.UIManager.getDefaults().getColor("Button.disabledToolBarBorderBackground"));
        jLabel18.setText("Quality Control Filters:");
        jLabel18.setOpaque(true);
        panelFiltersTitle.add(jLabel18, java.awt.BorderLayout.PAGE_START);

        panelOptions.add(panelFiltersTitle);
        panelFiltersTitle.setBounds(386, 13, 390, 265);

        panelInputs.add(panelOptions, java.awt.BorderLayout.CENTER);

        add(panelInputs, java.awt.BorderLayout.CENTER);
    }// </editor-fold>//GEN-END:initComponents

    private void selOutputDirButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_selOutputDirButtonActionPerformed
        if (evt.getSource() == selOutputDirButton) {
            fc.setCurrentDirectory(getLastDir());
            fc.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
            int returnVal = fc.showOpenDialog(ViewInputs.this);

            if (returnVal == JFileChooser.APPROVE_OPTION) {
                File file = fc.getSelectedFile();
                prefs.put(LAST_USED_FOLDER, fc.getSelectedFile().getAbsolutePath());
                //This is where a real application would open the file.
                fieldOutputDir.setText(file.getAbsolutePath());
                setEnabledInputs(true);
            }
        }
    }//GEN-LAST:event_selOutputDirButtonActionPerformed

    private void selPhenotypeButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_selPhenotypeButtonActionPerformed
        if (evt.getSource() == selPhenotypeButton) {
            fc.setCurrentDirectory(getLastDir());
            int returnVal = fc.showOpenDialog(ViewInputs.this);
            fc.setFileSelectionMode(JFileChooser.FILES_AND_DIRECTORIES);
            if (returnVal == JFileChooser.APPROVE_OPTION) {
                File file = fc.getSelectedFile();
                //This is where a real application would open the file.
                fieldPheno.setText(file.getAbsolutePath());
            }
        }
    }//GEN-LAST:event_selPhenotypeButtonActionPerformed

    private void selGenotypeButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_selGenotypeButtonActionPerformed
        //Handle open button action.
        if (evt.getSource() == selGenotypeButton) {
            fc.setCurrentDirectory(getLastDir());
            fc.setFileSelectionMode(JFileChooser.FILES_AND_DIRECTORIES);
            int returnVal = fc.showOpenDialog(ViewInputs.this);

            if (returnVal == JFileChooser.APPROVE_OPTION) {
                File file = fc.getSelectedFile();
                //This is where a real application would open the file.
                fieldGeno.setText(file.getAbsolutePath());
                prefs.put(LAST_USED_FOLDER, fc.getSelectedFile().getParent());
            }
        }        // TODO add your handling code here:
    }//GEN-LAST:event_selGenotypeButtonActionPerformed

    private void mapSelButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_mapSelButtonActionPerformed
         if (evt.getSource() == mapSelButton) {
            fc.setCurrentDirectory(getLastDir());
            int returnVal = fc.showOpenDialog(ViewInputs.this);
            fc.setFileSelectionMode(JFileChooser.FILES_AND_DIRECTORIES);
            if (returnVal == JFileChooser.APPROVE_OPTION) {
                File file = fc.getSelectedFile();
                //This is where a real application would open the file.
                fieldMap.setText(file.getAbsolutePath());
            }
        }
    }//GEN-LAST:event_mapSelButtonActionPerformed

    private void genotypeFormatCBoxActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_genotypeFormatCBoxActionPerformed
        String selection = genotypeFormatCBox.getSelectedItem().toString();
        System.out.println (selection);
        controller.onGenotypeFormat(selection);
    }//GEN-LAST:event_genotypeFormatCBoxActionPerformed

    private void fieldFilteringActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_fieldFilteringActionPerformed
        boolean flag = false;
        System.out.println(fieldFiltering.getSelectedItem().toString());
        if (fieldFiltering.getSelectedItem().toString().equals("FALSE")) {
            flag = false;
        } else {
            flag = true;
        }

        for (Component c : panelFilters.getComponents())
        c.setEnabled(flag);
    }//GEN-LAST:event_fieldFilteringActionPerformed

    private void fieldCorrectionActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_fieldCorrectionActionPerformed
        // TODO add your handling code here:
    }//GEN-LAST:event_fieldCorrectionActionPerformed

    private void fieldFilterMAFActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_fieldFilterMAFActionPerformed
        // TODO add your handling code here:
    }//GEN-LAST:event_fieldFilterMAFActionPerformed


    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JComboBox<String> fieldCorrection;
    private javax.swing.JTextField fieldFilterGENO;
    private javax.swing.JTextField fieldFilterHWE;
    private javax.swing.JTextField fieldFilterMAF;
    private javax.swing.JTextField fieldFilterMIND;
    private javax.swing.JComboBox<String> fieldFiltering;
    private javax.swing.JTextField fieldGeno;
    private javax.swing.JTextField fieldMap;
    private javax.swing.JComboBox<String> fieldModel;
    private javax.swing.JTextField fieldOutputDir;
    private javax.swing.JTextField fieldPheno;
    private javax.swing.JComboBox<String> fieldPloidy;
    private javax.swing.JTextField fieldR2;
    private javax.swing.JComboBox<String> fieldSNPs;
    private javax.swing.JTextField fieldSignificance;
    private javax.swing.Box.Filler filler1;
    private javax.swing.JComboBox<String> genotypeFormatCBox;
    private javax.swing.JLabel genotypeFormatLabel;
    private javax.swing.JLabel genotypeLabel;
    private javax.swing.JCheckBox jCheckBox1;
    private javax.swing.JLabel jLabel1;
    private javax.swing.JLabel jLabel15;
    private javax.swing.JLabel jLabel16;
    private javax.swing.JLabel jLabel17;
    private javax.swing.JLabel jLabel18;
    private javax.swing.JLabel jLabel2;
    private javax.swing.JLabel jLabel3;
    private javax.swing.JLabel jLabel4;
    private javax.swing.JLabel jLabel5;
    private javax.swing.JLabel jLabel6;
    private javax.swing.JLabel jLabel7;
    private javax.swing.JLabel jLabel9;
    private javax.swing.JTextField jTextField1;
    private javax.swing.JLabel labelGwasModel;
    private javax.swing.JLabel labelInputOutput;
    private javax.swing.JLabel labelOutputDir;
    private javax.swing.JLabel labelPloidy;
    private javax.swing.JLabel labelR2;
    private javax.swing.JLabel mapLabel;
    private javax.swing.JButton mapSelButton;
    private javax.swing.JPanel panelFilesTitle;
    private javax.swing.JPanel panelFilters;
    private javax.swing.JPanel panelFiltersTitle;
    private javax.swing.JInternalFrame panelGenotypeInfo;
    private javax.swing.JPanel panelInputs;
    private javax.swing.JPanel panelOptions;
    private javax.swing.JPanel panelParameters;
    private javax.swing.JPanel panelParametersTitle;
    private javax.swing.JPanel panelPaths;
    private javax.swing.JButton selGenotypeButton;
    private javax.swing.JButton selOutputDirButton;
    private javax.swing.JButton selPhenotypeButton;
    // End of variables declaration//GEN-END:variables
}
