package plugins;//package ij.plugin;

import ij.*;
import ij.gui.*;
import ij.measure.ResultsTable;
import ij.plugin.PlugIn;
import ij.plugin.filter.Analyzer;
import ij.plugin.filter.PlugInFilter;
import ij.process.*;
import java.text.NumberFormat;

/*
This plugin takes the Red and Green channel values from an RGB image or stack of images of cells, determines the
autophagy activity on a cell by cell basis for the give imagej outputing the results as a table with as many rows as
there are cells. For each cell the plugin computes :
    a)  whether or not there is colocalization for a given pixel, and stores the result in the Blue channel
        of the image or slice.  This method of storage for the colocalization data allows for the easy
        appreciation of the presence of colocalized pixels in a sample by superimposing the colocalization
        data without corrupting or modifying the original data.  For the colocalised pixels it measures the Integrated
        optical density.
    b)  area of the colocalised pixels expressed as #pixels.
    c)  area of the only RED pixels expressed as #pixels.
    d)  percentage of the colocalized pixels expressed as #Yellow pixels/#Total number of pixels.
    e)  ratio of red/green colocalized pixels expressed as #Yellow pixels/#Red pixels.

To alter the sensitivity of the selection algorithm, the user can specify the minimum ratio for the pixels in question,
as well as the threshold for the red and green channels.  Alternatively, the user can choose to have the threshold
values determined automatically (for each slice if a stack) via the getAutoThreshold() method in the ImageProcessor Class.
The colocalization data can be expressed as the average, max, or min of the corresponding red and green pixels, or as a
saturated pixel (255).  The colocalized data can also be displayed in a separate window as an 8bit or RGB image.  If any
data was present in the Blue channel of the original image, it will be replaced by the colocalization data.

					Elton Rexhepaj
					elton.rexhepaj@gmail.com
*/

public class AutophagyColocalisation implements PlugIn {

    private ImagePlus imp, imp2, imp3;
    private ImageProcessor ip;
    private ByteProcessor bp;
    private ColorProcessor cp,cp2,cp3;
    private ImageStack stack,stack2,stack3;
    private byte rpix[],gpix[],bpix[];
    private double ratio = 0.0;
    private double redthresh = 0;
    private double greenthresh = 0;
    private boolean atebit = true;
    private boolean auto = true;
    private int choice = 0;
    private int display = 2;
    private NumberFormat nf = NumberFormat.getInstance();
    private double prognum;
    private String progstr;
    private String otitle;
    private String titlestring;
    private ImageProcessor rip,gip;

    // 1 - add percentage of colocalised pixels
    private double percentColocPixels = 0;
    // 2 - add area of colocalised pixels
    private double areaColocPixels = 0;
    // 3 - add area of Red pixels
    private double areaRedPixels = 0;
    // 3 - add Integrated OD of colocalised pixels
    private double odColocPixels = 0;
    // 4 - add Ratio area (R+G)/R pixels
    private double ratioareaColocPixels1 = 0;
    // 5 - add Ratio area (R+G)/Total pixels
    private double ratioareaColocPixels2 = 0;

    public double getPercentColocPixels() {
        return percentColocPixels;
    }

    public void run(String arg){
        nf.setMaximumFractionDigits(2);
        nf.setMinimumFractionDigits(2);
        String errormessage;
        String choices[] = {"the average of the red and green","the max of the red and green",
                "the min of the red and green","a saturated pixel (255)"};
        String colocdisplay[] = {"Do not display",
                "an 8bit image","as an RGB image"};

        // create plugin Interface
        GenericDialog d = new GenericDialog("RG2B Colocalization", IJ.getInstance());
        //d.addCheckbox("Use auto thresholding", auto);
        d.addCheckbox("Use ImageJ to segment cells (DEFAULT=YES, Alternatively use CellProfiler)",auto);
        d.addMessage("or Manually Specify the following:");
        d.addNumericField("Minimum ratio between channels (0.0 to 1.0):", ratio, 2);
        d.addNumericField("Red channel lower threshold (0 to 255):", redthresh, 0);
        d.addNumericField("Green channel lower threshold (0 to 255):", greenthresh, 0);
        d.addMessage("Specify the output format:");
        d.addChoice("Set the colocalized pixels to:",choices,choices[1]);
        d.addChoice("Display the colocalization data separately as:",colocdisplay,colocdisplay[2]);
        d.showDialog();

        if(d.wasCanceled())return;
        auto = d.getNextBoolean();
        ratio = d.getNextNumber();
        if(auto)ratio=0;
        if(d.invalidNumber()) {
            errormessage = new String("Invalid value for ratio.\nAutomatic values will be calculated.");
            ratio = 0;
            IJ.showMessage("Error", errormessage);
        }
        redthresh = d.getNextNumber();
        if(d.invalidNumber()) {
            errormessage = new String("Invalid value for red channel threshold.\nAutomatic values will be calculated.");
            auto=true;
            ratio=0;
            IJ.showMessage("Error", errormessage);
        }
        greenthresh = d.getNextNumber();
        if(d.invalidNumber()) {
            errormessage = new String("Invalid value for green channel threshold.\nAutomatic values will be calculated.");
            auto=true;
            ratio=0;
            IJ.showMessage("Error", errormessage);
        }

        choice = d.getNextChoiceIndex();
        display = d.getNextChoiceIndex();
        if(ratio<0.0||ratio>1.0){
            ratio = 0;
            errormessage = new String("Invalid value for ratio.\nAutomatic values will be calculated.");
            IJ.showMessage(errormessage);
        }
        if(redthresh>255||redthresh<0){
            auto=true;
            errormessage = new String("Invalid value for red channel threshold.\nAutomatic values will be calculated.");
            IJ.showMessage(errormessage);
        }
        if(greenthresh>255||greenthresh<0){
            auto=true;
            errormessage = new String("Invalid value for green channel threshold.\nAutomatic values will be calculated.");
            IJ.showMessage(errormessage);
        }
        imp = IJ.getImage();
        otitle = imp.getTitle();
        if (imp.getBitDepth()!=24){
            IJ.showMessage("RGB Image required");
            return;
        }
        imp2 = new ImagePlus("RG2B Colocalization", imp.getStack());
        stack = imp.getStack();
        int m = imp.getWidth();
        int n = imp.getHeight();
        stack2 = imp2.createEmptyStack();
        if(display==1||display==2){
            imp3 = new ImagePlus("Colocalization Data", imp.getStack());
            stack3 = new ImageStack(m,n);
        }
        int dimension = m*n;
        int stacksize = stack.getSize();

        cp = (ColorProcessor)imp.getProcessor();
        rpix = new byte[dimension];
        gpix = new byte[dimension];
        bpix = new byte[dimension];

        for(int i=1;i<=stacksize;++i){
            cp = (ColorProcessor)stack.getProcessor(i);
            cp.getRGB(rpix,gpix,bpix);
            if(auto){
                rip = new ByteProcessor(m,n);
                rip.setPixels(rpix);
                redthresh = rip.getAutoThreshold();
                gip = new ByteProcessor(m,n);
                gip.setPixels(gpix);
                greenthresh = gip.getAutoThreshold();
            }

            // Run the colocalisation method to find colocalised pixels
            detectRGcolocalization(rpix,gpix,bpix);

            // create a new image to store the result of the colocalisation
            cp2 = new ColorProcessor(m,n);
            cp2.setRGB(rpix,gpix,bpix);
            stack2.addSlice(String.valueOf(i), cp2);

            if(display==1||display==2){
                cp3 = new ColorProcessor(m,n);
                cp3.setRGB(bpix,bpix,bpix);
                if(display==1){
                    ip = (new TypeConverter(cp3,true)).convertToByte();
                    stack3.addSlice(String.valueOf(i), ip);
                }
                else stack3.addSlice(String.valueOf(i), cp3);
            }
            prognum = (double)i/stacksize;
            progstr = nf.format(100*prognum);
            IJ.showStatus("RG2B Colocalization "+progstr+"% Done");
        }
        imp2.setStack("RG2B Colocalization - "+otitle, stack2);
        imp2.show();
        if(display==1||display==2){
            imp3.setStack("Colocalization Data - "+otitle, stack3);
            imp3.show();
        }

        // add results of the colocalisation plugin to the results table
        ResultsTable rt = Analyzer.getResultsTable();

        // create results table interface if interface not visible
        if (rt == null) {
            rt = new ResultsTable();
            Analyzer.setResultsTable(rt);
        }

        //TODO: Implement the code to retrive the cell masks
        int numbCells = 1;

        for (int i = 1; i <= numbCells; i++) {
            rt.incrementCounter();

            // 1 - add percentage of colocalised pixels
            rt.addValue("% RG Colocalised pixels", this.getPercentColocPixels());
            // 2 - add area of colocalised pixels
            rt.addValue("Area Colocalised pixels", this.getAreaColocPixels());
            // 3 - add Integrated OD of colocalised pixels
            rt.addValue("Area Red pixels", this.getAreaColocPixels());
            // 4 - add Integrated OD of colocalised pixels
            rt.addValue("Integrated OD pixels", this.getOdColocPixels());
            // 5 - add Ratio area (R+G)/R pixels
            rt.addValue("Ratio area (R+G)/R pixels", this.getRatioareaColocPixels1());
            // 6 - add Ratio area (R+G)/Total pixels
            rt.addValue("Ratio area (R+G)/Total pixels", this.getRatioareaColocPixels2());
        }

        // show(activate) the results window interface
        rt.show("Results");

    }


    private void detectRGcolocalization(byte redarray[], byte greenarray[], byte bluearray[]){
        int red, green;
        double ratio1, ratio2;
        for(int i=0;i<redarray.length;++i){
            red = 0xff & redarray[i];
            green = 0xff & greenarray[i];
            ratio1 = (double)green/(double)red;
            ratio2 = 1.0/ratio1;
            if(red>redthresh&&green>greenthresh&&ratio1>ratio&&ratio2>ratio){
                this.setAreaColocPixels(this.getAreaColocPixels()+1);
                this.setOdColocPixels(((red+green)/2)+this.getOdColocPixels());

                switch(choice){
                    case 0: bluearray[i]=(byte)((red+green)/2); break;
                    case 1: if(redarray[i]>greenarray[i])bluearray[i]=redarray[i];else bluearray[i]=greenarray[i];break;
                    case 2: if(redarray[i]<greenarray[i])bluearray[i]=redarray[i];else bluearray[i]=greenarray[i];break;
                    case 3: bluearray[i]=(byte)255;
                }
            }
            else {
                if(red>redthresh){
                    this.setAreaRedPixels(this.getAreaRedPixels()+1);
                }
                bluearray[i]=0;
            }
        }

        this.setPercentColocPixels(this.getAreaColocPixels()/(redarray.length));
        if(this.getAreaRedPixels()>0) {
            this.setRatioareaColocPixels1(this.getAreaColocPixels() / (this.getAreaRedPixels()));
        }
        if(this.getAreaRedPixels()>0) {
            this.setRatioareaColocPixels2(this.getAreaColocPixels() / (this.getAreaRedPixels()));
        }
    }

    public void setPercentColocPixels(double percentColocPixels) {
        this.percentColocPixels = percentColocPixels;
    }

    public double getAreaColocPixels() {
        return areaColocPixels;
    }

    public void setAreaColocPixels(double areaColocPixels) {
        this.areaColocPixels = areaColocPixels;
    }

    public double getAreaRedPixels() {
        return areaRedPixels;
    }

    public void setAreaRedPixels(double areaRedPixels) {
        this.areaRedPixels = areaRedPixels;
    }

    public double getOdColocPixels() {
        return odColocPixels;
    }

    public void setOdColocPixels(double odColocPixels) {
        this.odColocPixels = odColocPixels;
    }

    public double getRatioareaColocPixels1() {
        return ratioareaColocPixels1;
    }

    public void setRatioareaColocPixels1(double ratioareaColocPixels1) {
        this.ratioareaColocPixels1 = ratioareaColocPixels1;
    }

    public double getRatioareaColocPixels2() {
        return ratioareaColocPixels2;
    }

    public void setRatioareaColocPixels2(double ratioareaColocPixels2) {
        this.ratioareaColocPixels2 = ratioareaColocPixels2;
    }

}
