package jmultigwas;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Desktop;
import java.awt.Dimension;
import java.awt.Font;
import java.io.IOException;
import java.net.URISyntaxException;
import javax.swing.JEditorPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.event.HyperlinkEvent;
import javax.swing.event.HyperlinkListener;

public class ViewOutputs extends JPanel {
    Dimension dimension;
    JScrollPane scroll;
    JEditorPane outputArea;
    
    public ViewOutputs (Dimension d) {
        super ();
        setLayout (new BorderLayout ());
        dimension = new Dimension (d.width, (int) (0.9*d.height));       
    }

    public void init() {
        setBackground (Color.GRAY);
        // Init outputArea
        outputArea = new JEditorPane();
        outputArea.setEditable(false);//so its not editable
        outputArea.addHyperlinkListener(new HyperlinkListener() {
            public void hyperlinkUpdate(HyperlinkEvent e) {
                if(e.getEventType() == HyperlinkEvent.EventType.ACTIVATED) {
                    System.out.println(e.getURL());
                    // controller.browseFile (e.getURL);
                    
                    Desktop desktop = Desktop.getDesktop();
                            try {
                                desktop.browse(e.getURL().toURI());
                            } catch (IOException | URISyntaxException e1) {
                                e1.printStackTrace();
                            }
                }
            }
        });
        
        scroll = new JScrollPane(outputArea,
                JScrollPane.VERTICAL_SCROLLBAR_ALWAYS, JScrollPane.HORIZONTAL_SCROLLBAR_ALWAYS);
        scroll.setPreferredSize(dimension);
        add(scroll);
        outputArea.setFont(new Font("Monospaced",Font.PLAIN,15));
    }

    public void writeLine (String s, String type) {
        if (type=="html") {
            String text = outputArea.getText ();
            text = text.replaceAll("(\r\n|\n\r|\r|\n)", "<br />");
            text = text.replaceAll("(\t)", "&emsp;&emsp;&nbsp;&nbsp;");
            
            outputArea.setContentType("text/html");//set content as html
            outputArea.setText (text + "<br/><br/>"+ s);
        }else {
            outputArea.setText (outputArea.getText() +"\n"+s);
        }
        outputArea.setCaretPosition(outputArea.getDocument().getLength());
    }
    
    
    public void clearOutputs () {
        outputArea.setText ("");
        outputArea.setContentType("text/plain");//set content as html
    }
}
