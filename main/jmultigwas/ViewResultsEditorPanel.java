/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package jmultigwas;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.io.File;
import java.io.IOException;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.swing.JEditorPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;

public class ViewResultsEditorPanel extends JPanel {

    Dimension dimension;

    public ViewResultsEditorPanel(Dimension d) {
        super();
        setLayout(new BorderLayout());
        dimension = new Dimension(d.width, (int) (0.9 * d.height));
        this.setBackground(Color.LIGHT_GRAY);
        
    }
    
    public void showResults (String htmlFilename) {
        System.out.println (">>> Output report file: " + htmlFilename);
        try {
            File file = new File(htmlFilename);
            URL url = file.toURI().toURL();

            JEditorPane editorPane = new JEditorPane(url);
            editorPane.setContentType("text/html");
            //editorPane.setPage();
            editorPane.setEditable(false);

            JPanel noWrapPanel = new JPanel(new BorderLayout());
            noWrapPanel.add(editorPane);
            noWrapPanel.setPreferredSize(new Dimension(200, 200));

            //editorPane.setContentType("text/html");    
            JScrollPane scroll = new JScrollPane(noWrapPanel);
            scroll.setVerticalScrollBarPolicy(JScrollPane.VERTICAL_SCROLLBAR_ALWAYS);
            scroll.setViewportView(editorPane);
            scroll.setPreferredSize(dimension);
            add(scroll);
        } catch (MalformedURLException ex) {
            System.out.println(ex);
            Logger.getLogger(ViewResultsEditorPanel.class.getName()).log(Level.SEVERE, null, ex);
        } catch (IOException ex) {
            System.out.println(ex);
            Logger.getLogger(ViewResultsEditorPanel.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
}
