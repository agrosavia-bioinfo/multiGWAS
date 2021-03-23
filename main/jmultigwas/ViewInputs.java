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
        fieldPloidy.setSelectedIndex(0);
        genotypeText.setText("");
        fieldPheno.setText("");
        fieldSignificance.setText("0.05");
        fieldCorrection.setSelectedIndex(0);
        fieldModel.setSelectedIndex(0);
        fieldSNPs.setSelectedIndex(9);
        fieldFiltering.setSelectedIndex(0);
        fieldFilterMAF.setText("0.01");
        fieldFilterMIND.setText("0.1");
        fieldFilterGENO.setText("0.1");
        fieldFilterHWE.setText("1e-10");
    }

    public void clearInputs() {
        fieldPloidy.setSelectedIndex(0);
        fieldOutputDir.setText("");
        genotypeText.setText("");
        fieldPheno.setText("");
        mapText.setText("");
        fieldSignificance.setText("");
        fieldCorrection.setSelectedIndex(0);
        fieldModel.setSelectedIndex(0);
        fieldSNPs.setSelectedIndex(9);
        fieldFiltering.setSelectedIndex(0);
        fieldFilterMAF.setText("");
        fieldFilterMIND.setText("");
        fieldFilterGENO.setText("");
        fieldFilterHWE.setText("");
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
        selOutputDir.setEnabled(true);     
    }
    public void enableMapComponents (boolean enableFlag) {
        mapLabel.setEnabled (enableFlag);
        mapText.setEnabled (enableFlag);
        mapButton.setEnabled (enableFlag);
    } 
    
    // Init form with test mode options
    public void setTestMode (boolean testMode) {
        if (testMode) {
            genotypeText.setText ("/home/lg/agrosavia/GWAS-TOOL/development/tests/jmultiGWAS/example-genotype-tetra.csv");
            fieldPheno.setText ("/home/lg/agrosavia/GWAS-TOOL/development/tests/jmultiGWAS/example-phenotype.csv");
        }
    }

    public String getInputValues() {
        StringBuffer txt = new StringBuffer(500);
        String ln = System.lineSeparator();
        txt.append("default:" + ln);
        txt.append("  ploidy               : " + fieldPloidy.getSelectedItem().toString() + ln);
        txt.append("  genotypeFile         : " + '"' + genotypeText.getText() + '"' + ln);
        txt.append("  phenotypeFile        : " + '"' + fieldPheno.getText() + '"' + ln);
        txt.append("  genotypeFormat       : " + '"' + genotypeFormatCBox.getSelectedItem().toString() + '"' + ln);        
        txt.append("  mapFile              : " + '"' + mapText.getText() + '"' + ln); 
        txt.append("  significanceLevel    : " + fieldSignificance.getText() + ln);
        txt.append("  correctionMethod     : " + '"' + fieldCorrection.getSelectedItem().toString() + '"' + ln);
        txt.append("  gwasModel            : " + fieldModel.getSelectedItem().toString() + ln);
        txt.append("  nBest                : " + fieldSNPs.getSelectedItem().toString() + ln);
        txt.append("  filtering            : " + fieldFiltering.getSelectedItem().toString() + ln);
        txt.append("  MAF                  : " + fieldFilterMAF.getText() + ln);
        txt.append("  MIND                 : " + fieldFilterMIND.getText() + ln);
        txt.append("  GENO                 : " + fieldFilterGENO.getText() + ln);
        txt.append("  HWE                  : " + fieldFilterHWE.getText() + ln);
        txt.append("  tools                : " + '"' + controller.getToolsToRun() + '"' + ln);
        txt.append("  geneAction           : " + '"' + controller.getGeneAction() + '"' + ln);

        return txt.toString();
    }

    public boolean checkCompleteInfo() {
        if (fieldOutputDir.getText().equals ("")) return false;
        if (genotypeText.getText().equals ("")) return false;
        if (fieldPheno.getText().equals ("")) return false;
        if (fieldSignificance.getText().equals ("")) return false;
        if (fieldFilterMAF.getText().equals ("")) return false;
        if (fieldFilterMIND.getText().equals ("")) return false;
        if (fieldFilterGENO.getText().equals ("")) return false;
        if (fieldFilterHWE.getText().equals ("")) return false;
        
        if (genotypeFormatCBox.getSelectedItem().toString().equals ("KMatrix") ||
            genotypeFormatCBox.getSelectedItem().toString().equals ("FitPoly") ||
            genotypeFormatCBox.getSelectedItem().toString().equals ("Updog"))
                if (mapText.getText().equals("")) return false;
                
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
        genotypeText = new javax.swing.JTextField();
        genotypeSelButton = new javax.swing.JButton();
        genotypeFormatLabel = new javax.swing.JLabel();
        genotypeFormatCBox = new javax.swing.JComboBox<>();
        selOutputDir = new javax.swing.JButton();
        selPhenotypeBt = new javax.swing.JButton();
        mapLabel = new javax.swing.JLabel();
        mapText = new javax.swing.JTextField();
        mapButton = new javax.swing.JButton();
        panelOptions = new javax.swing.JPanel();
        panelParametersTitle = new javax.swing.JPanel();
        jLabel1 = new javax.swing.JLabel();
        panelParameters = new javax.swing.JPanel();
        jLabel8 = new javax.swing.JLabel();
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

        labelOutputDir.setText("Output Folder:");

        fieldOutputDir.setText("/home/lg/AAA");
        fieldOutputDir.setToolTipText("Select genotype file");

        fieldPheno.setToolTipText("Select phenotype file");
        fieldPheno.setPreferredSize(new java.awt.Dimension(90, 19));

        jLabel2.setText("Phenotype file:");

        genotypeLabel.setText("Genotype file:");

        genotypeText.setToolTipText("Select genotype file");
        genotypeText.setPreferredSize(new java.awt.Dimension(90, 19));

        genotypeSelButton.setText("Select...");
        genotypeSelButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                genotypeSelButtonActionPerformed(evt);
            }
        });

        genotypeFormatLabel.setText("Format:");

        genotypeFormatCBox.setModel(new javax.swing.DefaultComboBoxModel<>(new String[] { "GWASpoly", "KMatrix", "VCF", "FitPoly", "Updog", " " }));
        genotypeFormatCBox.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                genotypeFormatCBoxActionPerformed(evt);
            }
        });

        selOutputDir.setText("Select...");
        selOutputDir.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                selOutputDirActionPerformed(evt);
            }
        });

        selPhenotypeBt.setText("Select...");
        selPhenotypeBt.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                selPhenotypeBtActionPerformed(evt);
            }
        });

        mapLabel.setText("Map file:");
        mapLabel.setEnabled(false);

        mapText.setToolTipText("Select phenotype file");
        mapText.setEnabled(false);
        mapText.setPreferredSize(new java.awt.Dimension(90, 19));

        mapButton.setText("Select...");
        mapButton.setEnabled(false);
        mapButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                mapButtonActionPerformed(evt);
            }
        });

        javax.swing.GroupLayout panelPathsLayout = new javax.swing.GroupLayout(panelPaths);
        panelPaths.setLayout(panelPathsLayout);
        panelPathsLayout.setHorizontalGroup(
            panelPathsLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(panelPathsLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(panelPathsLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(labelOutputDir)
                    .addComponent(genotypeLabel))
                .addGap(35, 35, 35)
                .addGroup(panelPathsLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(fieldOutputDir)
                    .addComponent(genotypeText, javax.swing.GroupLayout.DEFAULT_SIZE, 297, Short.MAX_VALUE)
                    .addComponent(mapText, javax.swing.GroupLayout.DEFAULT_SIZE, 297, Short.MAX_VALUE)
                    .addComponent(fieldPheno, javax.swing.GroupLayout.DEFAULT_SIZE, 297, Short.MAX_VALUE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(panelPathsLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(selOutputDir, javax.swing.GroupLayout.Alignment.TRAILING)
                    .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, panelPathsLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                        .addComponent(selPhenotypeBt, javax.swing.GroupLayout.Alignment.TRAILING)
                        .addComponent(genotypeSelButton, javax.swing.GroupLayout.Alignment.TRAILING)
                        .addComponent(mapButton)))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(genotypeFormatLabel)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(genotypeFormatCBox, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addContainerGap())
            .addGroup(panelPathsLayout.createSequentialGroup()
                .addGap(12, 12, 12)
                .addGroup(panelPathsLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(mapLabel)
                    .addComponent(jLabel2))
                .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );
        panelPathsLayout.setVerticalGroup(
            panelPathsLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(panelPathsLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(panelPathsLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(selOutputDir, javax.swing.GroupLayout.PREFERRED_SIZE, 19, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(labelOutputDir)
                    .addComponent(fieldOutputDir, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addGroup(panelPathsLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(genotypeLabel)
                    .addComponent(genotypeText, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(genotypeSelButton, javax.swing.GroupLayout.PREFERRED_SIZE, 19, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(genotypeFormatLabel)
                    .addComponent(genotypeFormatCBox, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addGroup(panelPathsLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(mapLabel)
                    .addComponent(mapText, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(mapButton, javax.swing.GroupLayout.PREFERRED_SIZE, 19, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addGap(8, 8, 8)
                .addGroup(panelPathsLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(jLabel2)
                    .addComponent(fieldPheno, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(selPhenotypeBt, javax.swing.GroupLayout.PREFERRED_SIZE, 19, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );

        panelFilesTitle.add(panelPaths, java.awt.BorderLayout.CENTER);

        panelInputs.add(panelFilesTitle, java.awt.BorderLayout.NORTH);

        panelOptions.setBorder(new javax.swing.border.MatteBorder(null));

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

        jLabel8.setText("GWAS Model:");

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

        javax.swing.GroupLayout panelParametersLayout = new javax.swing.GroupLayout(panelParameters);
        panelParameters.setLayout(panelParametersLayout);
        panelParametersLayout.setHorizontalGroup(
            panelParametersLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(panelParametersLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(panelParametersLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(labelPloidy)
                    .addComponent(jLabel8)
                    .addComponent(jLabel3)
                    .addComponent(jLabel4)
                    .addComponent(jLabel5)
                    .addComponent(jLabel6))
                .addGap(5, 5, 5)
                .addGroup(panelParametersLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.CENTER)
                    .addComponent(fieldPloidy, 0, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addComponent(fieldModel, 0, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addComponent(fieldSignificance)
                    .addComponent(fieldCorrection, 0, 105, Short.MAX_VALUE)
                    .addComponent(fieldSNPs, 0, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addComponent(fieldFiltering, 0, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                .addGap(72, 72, 72))
        );
        panelParametersLayout.setVerticalGroup(
            panelParametersLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, panelParametersLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(panelParametersLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(fieldPloidy)
                    .addComponent(labelPloidy))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addGroup(panelParametersLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(jLabel8)
                    .addComponent(fieldModel))
                .addGap(8, 8, 8)
                .addGroup(panelParametersLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(jLabel3)
                    .addComponent(fieldSignificance))
                .addGap(18, 18, 18)
                .addGroup(panelParametersLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(jLabel4)
                    .addComponent(fieldCorrection))
                .addGap(11, 11, 11)
                .addGroup(panelParametersLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(fieldSNPs)
                    .addComponent(jLabel5))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addGroup(panelParametersLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(jLabel6)
                    .addComponent(fieldFiltering))
                .addGap(34, 34, 34))
        );

        panelParametersTitle.add(panelParameters, java.awt.BorderLayout.CENTER);

        panelOptions.add(panelParametersTitle);
        panelOptions.add(filler1);

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
                .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
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
                .addGap(22, 22, 22))
        );

        jLabel7.getAccessibleContext().setAccessibleName("");

        panelFiltersTitle.add(panelFilters, java.awt.BorderLayout.CENTER);

        jLabel18.setBackground(javax.swing.UIManager.getDefaults().getColor("Button.disabledToolBarBorderBackground"));
        jLabel18.setText("Quality Control Filters:");
        jLabel18.setOpaque(true);
        panelFiltersTitle.add(jLabel18, java.awt.BorderLayout.PAGE_START);

        panelOptions.add(panelFiltersTitle);

        panelInputs.add(panelOptions, java.awt.BorderLayout.WEST);

        add(panelInputs, java.awt.BorderLayout.CENTER);
    }// </editor-fold>//GEN-END:initComponents

    private void selOutputDirActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_selOutputDirActionPerformed
        if (evt.getSource() == selOutputDir) {
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
    }//GEN-LAST:event_selOutputDirActionPerformed

    private void selPhenotypeBtActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_selPhenotypeBtActionPerformed
        if (evt.getSource() == selPhenotypeBt) {
            fc.setCurrentDirectory(getLastDir());
            int returnVal = fc.showOpenDialog(ViewInputs.this);
            fc.setFileSelectionMode(JFileChooser.FILES_AND_DIRECTORIES);
            if (returnVal == JFileChooser.APPROVE_OPTION) {
                File file = fc.getSelectedFile();
                //This is where a real application would open the file.
                fieldPheno.setText(file.getAbsolutePath());
            }
        }
    }//GEN-LAST:event_selPhenotypeBtActionPerformed

    private void genotypeSelButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_genotypeSelButtonActionPerformed
        //Handle open button action.
        if (evt.getSource() == genotypeSelButton) {
            fc.setCurrentDirectory(getLastDir());
            fc.setFileSelectionMode(JFileChooser.FILES_AND_DIRECTORIES);
            int returnVal = fc.showOpenDialog(ViewInputs.this);

            if (returnVal == JFileChooser.APPROVE_OPTION) {
                File file = fc.getSelectedFile();
                //This is where a real application would open the file.
                genotypeText.setText(file.getAbsolutePath());
                prefs.put(LAST_USED_FOLDER, fc.getSelectedFile().getParent());
            }
        }        // TODO add your handling code here:
    }//GEN-LAST:event_genotypeSelButtonActionPerformed

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

    private void mapButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_mapButtonActionPerformed
        // TODO add your handling code here:
    }//GEN-LAST:event_mapButtonActionPerformed

    private void genotypeFormatCBoxActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_genotypeFormatCBoxActionPerformed
        String selection = genotypeFormatCBox.getSelectedItem().toString();
        System.out.println (selection);
        controller.onGenotypeFormat(selection);
    }//GEN-LAST:event_genotypeFormatCBoxActionPerformed

    private void fieldCorrectionActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_fieldCorrectionActionPerformed
        // TODO add your handling code here:
    }//GEN-LAST:event_fieldCorrectionActionPerformed


    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JComboBox<String> fieldCorrection;
    private javax.swing.JTextField fieldFilterGENO;
    private javax.swing.JTextField fieldFilterHWE;
    private javax.swing.JTextField fieldFilterMAF;
    private javax.swing.JTextField fieldFilterMIND;
    private javax.swing.JComboBox<String> fieldFiltering;
    private javax.swing.JComboBox<String> fieldModel;
    private javax.swing.JTextField fieldOutputDir;
    private javax.swing.JTextField fieldPheno;
    private javax.swing.JComboBox<String> fieldPloidy;
    private javax.swing.JComboBox<String> fieldSNPs;
    private javax.swing.JTextField fieldSignificance;
    private javax.swing.Box.Filler filler1;
    private javax.swing.JComboBox<String> genotypeFormatCBox;
    private javax.swing.JLabel genotypeFormatLabel;
    private javax.swing.JLabel genotypeLabel;
    private javax.swing.JButton genotypeSelButton;
    private javax.swing.JTextField genotypeText;
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
    private javax.swing.JLabel jLabel8;
    private javax.swing.JLabel labelInputOutput;
    private javax.swing.JLabel labelOutputDir;
    private javax.swing.JLabel labelPloidy;
    private javax.swing.JButton mapButton;
    private javax.swing.JLabel mapLabel;
    private javax.swing.JTextField mapText;
    private javax.swing.JPanel panelFilesTitle;
    private javax.swing.JPanel panelFilters;
    private javax.swing.JPanel panelFiltersTitle;
    private javax.swing.JPanel panelInputs;
    private javax.swing.JPanel panelOptions;
    private javax.swing.JPanel panelParameters;
    private javax.swing.JPanel panelParametersTitle;
    private javax.swing.JPanel panelPaths;
    private javax.swing.JButton selOutputDir;
    private javax.swing.JButton selPhenotypeBt;
    // End of variables declaration//GEN-END:variables
}
