package jmultigwas;

import java.awt.BorderLayout;
import java.awt.Desktop;
import java.io.File;
import java.io.IOException;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.swing.JFileChooser;
import javax.swing.JPanel;

public class ViewFiles extends JPanel {

    // Attributes
    Controller controller;
    
    JFileChooser fileChooser;

    public ViewFiles(Controller controller) {
        super();
        fileChooser = new JFileChooser(); 
        this.controller = controller;
        
        setLayout(new BorderLayout());
        add(fileChooser);
        
        
        fileChooser.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                fileChooser = (JFileChooser) evt.getSource();
                String command = evt.getActionCommand();
                if (command.equals(JFileChooser.APPROVE_SELECTION)) {
                    try {
                        File fileToOpen = fileChooser.getSelectedFile();
                        Desktop.getDesktop().open(fileToOpen);
                    } catch (IOException ex) {
                        Logger.getLogger(ViewFiles.class.getName()).log(Level.SEVERE, null, ex);
                    }
                } else if (command.equals(JFileChooser.CANCEL_SELECTION)) {
                    controller.onCancelButton ();
                }
            }
        });
    }
    public void changeDir (String dir) {
        System.out.println ("Changing to dir: " + dir);
        fileChooser.setCurrentDirectory(new File (dir));
    }
}
