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
import java.net.MalformedURLException;
import java.net.URL;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.swing.JPanel;

public class ViewResults extends JPanel {

    Dimension dimension;

    public ViewResults(Dimension d) {
        super();
        setLayout(new BorderLayout());
        dimension = new Dimension(d.width, (int) (0.9 * d.height));
        //setBackground(new java.awt.Color(153, 153, 255));
    }

    public void showResults(String htmlFilename) {
        try {
            File file = new File(htmlFilename);
            URL url = file.toURI().toURL();
            SwingFXWebView htmlView = new SwingFXWebView(url.toString(), dimension);

            add(htmlView);
        } catch (MalformedURLException ex) {
            System.out.println(ex);
            Logger.getLogger(ViewResults.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
}
