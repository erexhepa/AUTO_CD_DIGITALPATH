package ij.plugin;

/**
 * Created by rexhepaj on 01/12/15.
 */

import ij.*;
import ij.gui.GenericDialog;
import ij.gui.Overlay;
import ij.gui.Roi;
import ij.io.FileInfo;
import ij.io.OpenDialog;
import ij.io.Opener;
import ij.measure.Calibration;
import ij.measure.Measurements;
import ij.plugin.ContrastEnhancer;
import ij.plugin.FileInfoVirtualStack;
import ij.plugin.FolderOpener;
import ij.process.ImageProcessor;
import ij.process.ImageStatistics;
import ij.util.StringSorter;
import ij.process.StackStatistics;
import ij.plugin.PlugIn;
import java.awt.List;
import java.lang.String;
import java.awt.*;
import java.awt.event.ItemEvent;
import java.awt.event.TextEvent;
import java.awt.image.ColorModel;
import java.io.File;
import ij.process.FloatProcessor;
import ij.process.FloatStatistics;

import java.io.*;
import java.util.*;
import ij.io.FileSaver;

/** Implements the File/Import/Image Sequence command, which
 opens a folder of images as a stack. */
public class BFX_TIFconverter implements PlugIn,Measurements {

    private static String[] excludedTypes = {".txt", ".lut", ".roi", ".pty", ".hdr", ".java", ".ijm", ".py", ".js", ".bsh", ".xml"};
    private static boolean staticSortFileNames = true;
    private static boolean staticOpenAsVirtualStack;
    private boolean convertToRGB = false;
    private boolean sortFileNames = true;
    private boolean openAsVirtualStack = true;
    private double scale = 100.0;
    private int n, start, increment;
    private String filter;
    private String legacyRegex;
    private FileInfo fi;
    private String info1;
    private ImagePlus image;
    private boolean saveImage = false;
    private long t0;
    private ImageStack[] calibrationSet;
    private double[][] scalingParam ;
    private String outdir = "/ac-ulm/transfert/Plateforme Biophenics/plateformes/incell/RUBIES_JPEG";
    private String indir  = "";
    private List dirList ;
    private Set<String> listChannelsUnique ;


    /** Opens the images in the specified directory as a stack. Displays
     directory chooser and options dialogs if the argument is null. */
    public static ImagePlus open(String path) {
        BFX_TIFconverter fo = new BFX_TIFconverter();
        fo.saveImage = true;
        fo.run(path);
        return fo.image;
    }

    /** Opens the images in the specified directory as a stack. Displays
     directory chooser and options dialogs if the argument is null. */
    public ImagePlus openFolder(String path) {
        saveImage = true;
        run(path);
        return image;
    }

    public String[] trimStrings(String[] listToTrim) {
        String[] trimmedList = new String[listToTrim.length];

        int listLen = listToTrim.length;

        for ( int i = 0; i < listLen; i++ ) {
            String currentFname = listToTrim[i].trim();
            int strLen = currentFname.length();
            // remove changing bits from the filename and add to the list
            trimmedList[i] = currentFname.substring(currentFname.length() - 8, currentFname.length());
        }

        return trimmedList;
    }

    public ImageStack loadImsequenceFromList(String[] list, String directory){
        String title = directory;

        if (title.endsWith(File.separator) || title.endsWith("/"))
            title = title.substring(0, title.length()-1);
        int index = title.lastIndexOf(File.separatorChar);
        if (index!=-1) title = title.substring(index + 1);
        if (title.endsWith(":"))
            title = title.substring(0, title.length()-1);

        IJ.register(FolderOpener.class);
        list = trimFileList(list);
        if (list==null) return null;
        if (IJ.debugMode) IJ.log("FolderOpener: "+directory+" ("+list.length+" files)");
        int width=0, height=0, stackSize=1, bitDepth=0;
        ImageStack stack = null;
        double min = Double.MAX_VALUE;
        double max = -Double.MAX_VALUE;
        Calibration cal = null;
        boolean allSameCalibration = true;
        IJ.resetEscape();
        Overlay overlay = null;
        n = list.length;
        start = 1;
        increment = 1;
        try {

            for (int i=0; i<list.length; i++) {
                Opener opener = new Opener();
                opener.setSilentMode(true);
                IJ.redirectErrorMessages(true);
                ImagePlus imp = opener.openImage(directory, list[i]);
                IJ.redirectErrorMessages(false);
                if (imp!=null) {
                    width = imp.getWidth();
                    height = imp.getHeight();
                    bitDepth = imp.getBitDepth();

                    break;
                }
            }//
            if (width==0) {
                IJ.error("Sequence Reader", "This folder does not appear to contain\n"
                        + "any TIFF, JPEG, BMP, DICOM, GIF, FITS or PGM files.\n \n"
                        + "   \""+directory+"\"");
                return null;
            }


            String pluginName = "Sequence Reader";
            if (legacyRegex!=null)
                pluginName += "(legacy)";
            list = getFilteredList(list, filter, pluginName);
            if (list==null)
                return null;
            IJ.showStatus("");
            t0 = System.currentTimeMillis();
            if (sortFileNames)
                list = StringSorter.sortNumerically(list);

            if (n<1)
                n = list.length;
            if (start<1 || start>list.length)
                start = 1;
            if (start+n-1>list.length)
                n = list.length-start+1;
            int count = 0;
            int counter = 0;
            ImagePlus imp = null;
            boolean firstMessage = true;
            boolean fileInfoStack = false;
            for (int i=start-1; i<list.length; i++) {
                if ((counter++%increment)!=0)
                    continue;
                Opener opener = new Opener();
                opener.setSilentMode(true);
                IJ.redirectErrorMessages(true);
                if ("RoiSet.zip".equals(list[i])) {
                    IJ.open(directory+list[i]);
                    imp = null;
                } else if (!openAsVirtualStack||stack==null) {
                    imp = opener.openImage(directory, list[i]);
                    stackSize = imp!=null?imp.getStackSize():1;
                }
                IJ.redirectErrorMessages(false);
                if (imp!=null && stack==null) {
                    width = imp.getWidth();
                    height = imp.getHeight();
                    bitDepth = imp.getBitDepth();
                    fi = imp.getOriginalFileInfo();
                    ImageProcessor ip = imp.getProcessor();
                    min = ip.getMin();
                    max = ip.getMax();
                    cal = imp.getCalibration();
                    if (convertToRGB) bitDepth = 24;
                    ColorModel cm = imp.getProcessor().getColorModel();
                    if (openAsVirtualStack) {
                        if (stackSize>1) {
                            stack = new FileInfoVirtualStack();
                            fileInfoStack = true;
                        } else
                            stack = new VirtualStack(width, height, cm, directory);
                        ((VirtualStack)stack).setBitDepth(bitDepth);
                    } else if (scale<100.0)
                        stack = new ImageStack((int)(width*scale/100.0), (int)(height*scale/100.0), cm);
                    else
                        stack = new ImageStack(width, height, cm);
                    info1 = (String)imp.getProperty("Info");
                }
                if (imp==null)
                    continue;
                if (imp.getWidth()!=width || imp.getHeight()!=height) {
                    IJ.log(list[i] + ": wrong size; "+width+"x"+height+" expected, "+imp.getWidth()+"x"+imp.getHeight()+" found");
                    continue;
                }
                String label = imp.getTitle();
                if (stackSize==1) {
                    String info = (String)imp.getProperty("Info");
                    if (info!=null)
                        label += "\n" + info;
                }
                if (imp.getCalibration().pixelWidth!=cal.pixelWidth)
                    allSameCalibration = false;
                ImageStack inputStack = imp.getStack();
                Overlay overlay2 = imp.getOverlay();
                if (overlay2!=null && !openAsVirtualStack) {
                    if (overlay==null)
                        overlay = new Overlay();
                    for (int j=0; j<overlay2.size(); j++) {
                        Roi roi = overlay2.get(j);
                        int position = roi.getPosition();
                        if (position==0)
                            roi.setPosition(count+1);
                        overlay.add(roi);
                    }
                }

                if (openAsVirtualStack) {
                    if (fileInfoStack)
                        openAsFileInfoStack((FileInfoVirtualStack)stack, directory+list[i]);
                    else
                        ((VirtualStack)stack).addSlice(list[i]);
                } else {
                    for (int slice=1; slice<=stackSize; slice++) {
                        int bitDepth2 = imp.getBitDepth();
                        String label2 = label;
                        ImageProcessor ip = null;
                        if (stackSize>1) {
                            String sliceLabel = inputStack.getSliceLabel(slice);
                            if (sliceLabel!=null)
                                label2=sliceLabel;
                            else if (label2!=null && !label2.equals(""))
                                label2 += ":"+slice;
                        }
                        ip = inputStack.getProcessor(slice);
                        if (convertToRGB) {
                            ip = ip.convertToRGB();
                            bitDepth2 = 24;
                        }
                        if (bitDepth2!=bitDepth) {
                            if (bitDepth==8 && bitDepth2==24) {
                                ip = ip.convertToByte(true);
                                bitDepth2 = 8;
                            } else if (bitDepth==32) {
                                ip = ip.convertToFloat();
                                bitDepth2 = 32;
                            } else if (bitDepth==24) {
                                ip = ip.convertToRGB();
                                bitDepth2 = 24;
                            }
                        }
                        if (bitDepth2!=bitDepth) {
                            IJ.log(list[i] + ": wrong bit depth; "+bitDepth+" expected, "+bitDepth2+" found");
                            break;
                        }
                        if (scale<100.0)
                            ip = ip.resize((int)(width*scale/100.0), (int)(height*scale/100.0));
                        if (ip.getMin()<min) min = ip.getMin();
                        if (ip.getMax()>max) max = ip.getMax();
                        stack.addSlice(label2, ip);
                    }
                }
                count++;
                //IJ.showStatus(count+"/"+n);
                //IJ.showProgress(count, n);
                if (count>=n)
                    break;
                if (IJ.escapePressed())
                {IJ.beep(); break;}
            }
        } catch(OutOfMemoryError e) {
            IJ.outOfMemory("FolderOpener");
            if (stack!=null) stack.trim();
        }


        return stack;
    }


    public void walk( String path ) {

        File root = new File( path );
        File[] list = root.listFiles();

        if (list == null) return;

        for ( File f : list ) {
            if ( f.isDirectory() ) {
                walk( f.getAbsolutePath() );
                System.out.println( "Dir:" + f.getAbsoluteFile() );
            }
            else {
                System.out.println( "File:" + f.getAbsoluteFile() );
            }
        }
    }

    /**
     * Method to scan folder automatically and detect how many channels are present and the scaling factor for each
     * respective channel.
     */
    public void setConvertionParam(String directory){
        // get directory listing of all files ending with '.tif'
        String[] list = (new File(directory)).list(new FilenameFilter() {
            public boolean accept(File dir, String name) {
                return name.toLowerCase().endsWith(".tif");
            }
        });

        // create copy of directory listing with only the channel name
        String[] listChannels = this.trimStrings(list);

        // get unique channels
        java.util.List<String> wordList = Arrays.asList(listChannels);
        listChannelsUnique = new HashSet<String>( wordList);

        //int[][] scalingParam = new int [listChannelsUnique.size()][2];
        calibrationSet = new ImageStack[listChannelsUnique.size()];
        scalingParam   = new double [listChannelsUnique.size()][2];

        final Object[] arrayChannels = listChannelsUnique.toArray();

        // initialize to default values the scaling parameters
        for (int chanIndx=0;chanIndx<listChannelsUnique.size();chanIndx++){
            scalingParam[chanIndx][0] = 100;
            scalingParam[chanIndx][1] = 350;
            // get list of each channel type
            final String finalChanIndx = arrayChannels[chanIndx].toString();

            String[] listChannel = (new File(directory)).list(new FilenameFilter() {
                public boolean accept(File dir, String name) {
                    return name.endsWith(finalChanIndx);
                }
            });

            // get random selection of 10 elements

            Random myRandomizer = new Random();
            int randSeleCalibration = 10;
            String[] listChannelRand;

            if(listChannel.length<=10){
                randSeleCalibration = listChannel.length;
                listChannelRand = new String[listChannel.length];
            }else{
                listChannelRand = new String[10];
            }

            for(int randIndxSel=0;randIndxSel<randSeleCalibration;randIndxSel++){
                listChannelRand[randIndxSel] = listChannel[myRandomizer.nextInt(listChannel.length)];
            }

            calibrationSet[chanIndx] = this.loadImsequenceFromList(listChannelRand, directory);
            String stackTitle = "Processing calibration image set - ";
            stackTitle.concat(finalChanIndx);

            ImagePlus impProcess = new ImagePlus(stackTitle, calibrationSet[chanIndx]);
            ImageProcessor ip = impProcess.getProcessor();
            double minStack = ip.getMin();
            double maxStack = ip.getMax();

            StackStatistics sps = new StackStatistics(impProcess, 512, minStack, maxStack);
            double spsMin   = sps.min;
            double spsMax   = sps.max;
            int indxMaxpix = getPercentile(sps.getHistogram(), (float) (0.99));

            scalingParam[chanIndx][0] = sps.min;
            scalingParam[chanIndx][1] = sps.binSize*(indxMaxpix+2)+sps.min;

            ImageStack stack = impProcess.getStack();
            Calibration cal = impProcess.getCalibration();
            ImageStatistics stats = ImageStatistics.getStatistics(ip, MEAN, cal);

            //System.out.print("Debug");
            //
            //IJ.showProgress(chanIndx,listChannelsUnique.size());
        }

    }

    public int getPercentile(long[] imHist, float percThresh){
        int percentReturn = 0;
        long[] filteredHist = new long[imHist.length-2];
        float[] filteredHistpercent = new float[imHist.length-2];

        // remove first bin count (i.e. background)
        for(int i=0;i<filteredHist.length;i++){
            filteredHist[i] = imHist[i+2];
        }

        // normalise histCounts
        long normFactor = 0;

        for(int i=0;i<filteredHist.length;i++){
            normFactor = normFactor+filteredHist[i];
        }

        for(int i=0;i<filteredHist.length;i++){
            filteredHistpercent[i] = (float)filteredHist[i]/(float)normFactor;
        }

        float cmdSum = 0;
        for(int i=0;i<filteredHistpercent.length;i++){
            cmdSum = cmdSum + filteredHistpercent[i];
            if(cmdSum>=percThresh){
                percentReturn = i;
                break;
            }
        }

        return percentReturn;
    }

    public void convertBFXfolder(String directory, String outDirectory){

        final Object[] arrayChannels = listChannelsUnique.toArray();

        for (int chanIndx=0;chanIndx<listChannelsUnique.size();chanIndx++){
            //scalingParam[chanIndx][0] = 100;
            //scalingParam[chanIndx][1] = 350;

            // get list of each channel type
            final String finalChanIndx = arrayChannels[chanIndx].toString();

            String[] listChannel = (new File(directory)).list(new FilenameFilter() {
                public boolean accept(File dir, String name) {
                    return name.endsWith(finalChanIndx);
                }
            });

            // iterate through all the images of the specified channel
            for(int chanFileIndx=0;chanFileIndx<listChannel.length;chanFileIndx++){
                ImagePlus imp = null;
                ContrastEnhancer contrastEnhancer = new ContrastEnhancer();
                Opener opener = new Opener();
                opener.setSilentMode(false);
                IJ.redirectErrorMessages(true);

                imp = opener.openImage(directory, listChannel[chanFileIndx]);
                ImageProcessor ip = imp.getProcessor();
                double min = ip.getMin();
                double max = ip.getMax();
                Calibration cal = imp.getCalibration();
                ip.setMinAndMax(scalingParam[chanIndx][0],scalingParam[chanIndx][1]);
                imp.show();
                IJ.showProgress(chanIndx*listChannel.length+chanFileIndx,listChannel.length*listChannelsUnique.size());
                System.out.println("Converting File :: "+listChannel[chanFileIndx]);



                File theDir = new File(outDirectory);

                // if the directory does not exist, create it
                if (!theDir.exists()) {
                    System.out.println("creating directory: " + outDirectory);
                    boolean result = false;

                    try{
                        theDir.mkdir();
                        result = true;
                    }
                    catch(SecurityException se){
                        //handle it
                    }
                    if(result) {
                        System.out.println("DIR created");
                    }
                }


                //
                String imfOutpath = outDirectory + "/" +listChannel[chanFileIndx].replace(".tif",".jpeg");
                saveCustomBFX(imp,imfOutpath);
                imp.close();
            }
        }

    }

    public void run(String arg) {
        boolean isMacro = Macro.getOptions()!=null;
        String directory = null;
        if (arg!=null && !arg.equals("")) {
            directory = arg;
        } else {
            if (!isMacro) {
                sortFileNames = staticSortFileNames;
                openAsVirtualStack = staticOpenAsVirtualStack;
            }
            arg = null;
            String title = "Open Image Sequence...";
            String macroOptions = Macro.getOptions();
            if (macroOptions!=null) {
                directory = Macro.getValue(macroOptions, title, null);
                if (directory!=null) {
                    directory = OpenDialog.lookupPathVariable(directory);
                    File f = new File(directory);
                    if (!f.isDirectory() && (f.exists()||directory.lastIndexOf(".")>directory.length()-5))
                        directory = f.getParent();
                }
                legacyRegex = Macro.getValue(macroOptions, "or", "");
                if (legacyRegex.equals(""))
                    legacyRegex = null;
            }
            if (directory==null) {
                directory = IJ.getDirectory("Select input directory");
                indir     = directory;
                outdir    = IJ.getDirectory("Select output directory");
            }
        }
        if (directory==null)
            return;

        // get directory listing of all files ending with '.tif'
        String[] listDir = (new File(directory)).list(new FilenameFilter() {
            public boolean accept(File dir, String name) {
                return name.toLowerCase().matches(".*[0-9]{8,8}.*");
            }
        });

        // Iterate through all directories and complete scaling calibration and convertion
        for(int dirIndx=0;dirIndx<listDir.length;dirIndx++){
            // get scaling parameters for plate-i
            // TODO: Implement
            this.setConvertionParam(indir+"/"+listDir[dirIndx]);

            // convert all images for all channels in plate-i
            // TODO: Implement
            String barCodePlate = "12345678";
            barCodePlate = getBarcodeFromDir(listDir[dirIndx]);
            this.convertBFXfolder(indir+listDir[dirIndx],outdir+barCodePlate);
        }

        if (listDir==null)
            return;
    }

    public String getBarcodeFromDir(String dirFname){
        String barCode =  "12345678";

        // match a regular expression of 8 numerical digits sequence in the directory name

        if(dirFname.matches(".*[0-9]{8,8}.*")){
            String[] barCodeTrimSeq = dirFname.split("[0-9]{8,8}");

            if(barCodeTrimSeq.length>0){
                for(int trimIndxSeq=0;trimIndxSeq<barCodeTrimSeq.length;trimIndxSeq++){
                    dirFname = dirFname.replace(barCodeTrimSeq[trimIndxSeq],"");
                }
                barCode = dirFname;
            }else{
                barCode = dirFname;
            }
        }

        return barCode;
    }

    public void saveCustomBFX(ImagePlus imp, String path) {

        try {
            FileSaver fs = new FileSaver(imp);
            fs.saveAsJpeg(path);

        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    private void openAsFileInfoStack(FileInfoVirtualStack stack, String path) {
        FileInfo[] info = Opener.getTiffFileInfo(path);
        if (info==null || info.length==0)
            return;
        int n =info[0].nImages;
        if (info.length==1 && n>1) {
            long size = fi.width*fi.height*fi.getBytesPerPixel();
            for (int i=0; i<n; i++) {
                FileInfo fi = (FileInfo)info[0].clone();
                fi.nImages = 1;
                fi.longOffset = fi.getOffset() + i*(size + fi.gapBetweenImages);
                stack.addImage(fi);
            }
        } else
            stack.addImage(info[0]);
    }

    boolean showDialog(ImagePlus imp, String[] list) {
        int fileCount = list.length;
        FolderOpenerDialog gd = new FolderOpenerDialog("Sequence Options", imp, list);
        gd.addNumericField("Number of images:", fileCount, 0);
        gd.addNumericField("Starting image:", 1, 0);
        gd.addNumericField("Increment:", 1, 0);
        gd.addNumericField("Scale images:", scale, 0, 4, "%");
        gd.addStringField("File name contains:", "", 10);
        gd.setInsets(0,45,0);
        gd.addMessage("(enclose regex in parens)", null, Color.darkGray);
        gd.addCheckbox("Convert_to_RGB", convertToRGB);
        gd.addCheckbox("Sort names numerically", sortFileNames);
        gd.addCheckbox("Use virtual stack", openAsVirtualStack);
        gd.addMessage("10000 x 10000 x 1000 (100.3MB)");
        gd.addHelp(IJ.URL+"/docs/menus/file.html#seq1");
        gd.setSmartRecording(true);
        gd.showDialog();
        if (gd.wasCanceled())
            return false;
        n = (int)gd.getNextNumber();
        start = (int)gd.getNextNumber();
        increment = (int)gd.getNextNumber();
        if (increment<1)
            increment = 1;
        scale = gd.getNextNumber();
        if (scale<5.0) scale = 5.0;
        if (scale>100.0) scale = 100.0;
        filter = gd.getNextString();
        if (legacyRegex!=null)
            filter = "("+legacyRegex+")";
        convertToRGB = gd.getNextBoolean();
        sortFileNames = gd.getNextBoolean();
        openAsVirtualStack = gd.getNextBoolean();
        if (openAsVirtualStack)
            scale = 100.0;
        if (!IJ.macroRunning()) {
            staticSortFileNames = sortFileNames;
            staticOpenAsVirtualStack = openAsVirtualStack;
        }
        return true;
    }

    public static String[] getFilteredList(String[] list, String filter, String title) {
        boolean isRegex = false;
        if (filter!=null && (filter.equals("") || filter.equals("*")))
            filter = null;
        if (list==null || filter==null)
            return list;
        if (title==null) {
            String[] list2 = new String[list.length];
            for (int i=0; i<list.length; i++)
                list2[i] = list[i];
            list = list2;
        }

        if (filter.length()>=2 && filter.startsWith("(")&&filter.endsWith(")")) {
            filter = filter.substring(1,filter.length()-1);
            isRegex = true;
        }
        int filteredImages = 0;
        for (int i=0; i<list.length; i++) {
            if (isRegex && containsRegex(list[i],filter,title!=null&&title.contains("(legacy)")))
                filteredImages++;
            else if (list[i].contains(filter))
                filteredImages++;
            else
                list[i] = null;
        }
        if (filteredImages==0) {
            if (title!=null) {
                if (isRegex)
                    IJ.error(title, "None of the file names contain the regular expression '"+filter+"'.");
                else
                    IJ.error(title, "None of the "+list.length+" files contain '"+filter+"' in the name.");
            }
            return null;
        }
        String[] list2 = new String[filteredImages];
        int j = 0;
        for (int i=0; i<list.length; i++) {
            if (list[i]!=null)
                list2[j++] = list[i];
        }
        list = list2;
        return list;
    }

    private static boolean containsRegex(String name, String regex, boolean legacy) {
        boolean contains = false;
        try {
            if (legacy)
                contains = name.matches(regex);
            else
                contains = name.replaceAll(regex,"").length()!=name.length();
            IJ.showStatus("");
        } catch(Exception e) {
            String msg = e.getMessage();
            int index = msg.indexOf("\n");
            if (index>0)
                msg = msg.substring(0,index);
            IJ.showStatus("Regex error: "+msg);
            contains = true;
        }
        return contains;
    }

    /** Removes names that start with "." or end with ".db", ".txt", ".lut", "roi", ".pty", ".hdr", ".py", etc. */
    public String[] trimFileList(String[] rawlist) {
        int count = 0;
        for (int i=0; i< rawlist.length; i++) {
            String name = rawlist[i];
            if (name.startsWith(".")||name.equals("Thumbs.db")||excludedFileType(name))
                rawlist[i] = null;
            else
                count++;
        }
        if (count==0) return null;
        String[] list = rawlist;
        if (count<rawlist.length) {
            list = new String[count];
            int index = 0;
            for (int i=0; i< rawlist.length; i++) {
                if (rawlist[i]!=null)
                    list[index++] = rawlist[i];
            }
        }
        return list;
    }

    /* Returns true if 'name' ends with ".txt", ".lut", ".roi", ".pty", ".hdr", ".java", ".ijm", ".py", ".js" or ".bsh. */
    public static boolean excludedFileType(String name) {
        if (name==null) return true;
        for (int i=0; i<excludedTypes.length; i++) {
            if (name.endsWith(excludedTypes[i]))
                return true;
        }
        return false;
    }

    public void openAsVirtualStack(boolean b) {
        openAsVirtualStack = b;
    }

    public void sortFileNames(boolean b) {
        sortFileNames = b;
    }

    /** Sorts file names containing numerical components.
     * @see ij.util.StringSorter#sortNumerically
     * Author: Norbert Vischer
     */
    public String[] sortFileList(String[] list) {
        return StringSorter.sortNumerically(list);
    }

    class FolderOpenerDialog extends GenericDialog {
        ImagePlus imp;
        int fileCount;
        boolean eightBits, rgb;
        String[] list;
        //boolean isRegex;

        public FolderOpenerDialog(String title, ImagePlus imp, String[] list) {
            super(title);
            this.imp = imp;
            this.list = list;
            this.fileCount = list.length;
        }

        protected void setup() {
            eightBits = ((Checkbox)checkbox.elementAt(0)).getState();
            rgb = ((Checkbox)checkbox.elementAt(1)).getState();
            setStackInfo();
        }

        public void itemStateChanged(ItemEvent e) {
        }

        public void textValueChanged(TextEvent e) {
            setStackInfo();
        }

        void setStackInfo() {
            if (imp==null)
                return;
            int width = imp.getWidth();
            int height = imp.getHeight();
            int depth = imp.getStackSize();
            int bytesPerPixel = 1;
            int n = getNumber(numberField.elementAt(0));
            int start = getNumber(numberField.elementAt(1));
            int inc = getNumber(numberField.elementAt(2));
            double scale = getNumber(numberField.elementAt(3));
            if (scale<5.0) scale = 5.0;
            if (scale>100.0) scale = 100.0;
            if (n<1) n = fileCount;
            if (start<1 || start>fileCount) start = 1;
            if (start+n-1>fileCount)
                n = fileCount-start+1;
            if (inc<1) inc = 1;
            TextField tf = (TextField)stringField.elementAt(0);
            String filter = tf.getText();
            int n3 = Integer.MAX_VALUE;
            String[] filteredList = getFilteredList(list, filter, null);
            if (filteredList!=null)
                n3 = filteredList.length;
            else
                n3 = 0;
            if (n3<n)
                n = n3;
            switch (imp.getType()) {
                case ImagePlus.GRAY16:
                    bytesPerPixel=2;break;
                case ImagePlus.COLOR_RGB:
                case ImagePlus.GRAY32:
                    bytesPerPixel=4; break;
            }
            if (eightBits)
                bytesPerPixel = 1;
            if (rgb)
                bytesPerPixel = 4;
            width = (int)(width*scale/100.0);
            height = (int)(height*scale/100.0);
            int n2 = ((fileCount-start+1)*depth)/inc;
            if (n2<0) n2 = 0;
            if (n2>n) n2 = n;
            double size = ((double)width*height*n2*bytesPerPixel)/(1024*1024);
            ((Label)theLabel).setText(width+" x "+height+" x "+n2+" ("+IJ.d2s(size,1)+"MB)");
        }

        public int getNumber(Object field) {
            TextField tf = (TextField)field;
            String theText = tf.getText();
            double value;
            Double d;
            try {d = new Double(theText);}
            catch (NumberFormatException e){
                d = null;
            }
            if (d!=null)
                return (int)d.doubleValue();
            else
                return 0;
        }

    } // FolderOpenerDialog

} // FolderOpener

