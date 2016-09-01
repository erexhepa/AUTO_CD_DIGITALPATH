package AUTO_COL_DECONV;

/***************************************************
 *
 * License
 *
 * PlugIn:		CDReview
 *
 * File: 		CDR_run.java
 *
 * Description:	CDReview CDR_run - Software for
 * 				testing and comparison of colour deconvolution methods
 * 				and for creation of benchmark RGB colour images
 *
 *				This software is provided as a supplement to:
 *
 *				Haub, P. and Meckel, T. A Model based Survey of Colour Deconvolution
 *				in Diagnostic Brightfield Microscopy: Error Estimation and Spectral Consideration.
 *				Sci. Rep. 5, 12096; doi: 10.1038/srep12096 (2015).
 *
 *
 * @author: 	Peter Haub, 2015
 *
 * Copyright(C) 2015 Peter Haub
 *              phaub@dipsystems.de
 *              www.dipsystems.de
 *
 * Version:     0.1.0 - July 2015
 *
 * License:		This program is free software; you can redistribute it and/or modify it
 * 				under the terms of the GNU General Public License as published by the
 * 				Free Software Foundation; either version 2 of the License, or (at your
 * 				option) any later version.
 *
 * 				This program is distributed in the hope that it will be useful, but
 * 				WITHOUT ANY WARRANTY; without even the implied warranty of
 * 				MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * 				General Public License for more details.
 *
 * 				You should have received a copy of the GNU General Public License along
 * 				with this program; if not, write to the Free Software Foundation, Inc.,
 * 				51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 *
 ***************************************************/


import ij.*;
import ij.gui.GenericDialog;
import ij.gui.NewImage;
import ij.plugin.PlugIn;
import ij.process.ColorProcessor;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import ij.process.StackConverter;
import org.apache.commons.math3.linear.*;

import java.io.File;
import java.util.ArrayList;


public class CDR_run_cmdline implements PlugIn {

    static int R = 0;
    static int G = 1;
    static int B = 2;

    String[] modes2Stain = { "GB", "GR", "RB"};
    String[] modesDecon = { "DET", "LU", "QR"};

    ImagePlus impRGB, impDecon;
    double[] aR, aG, aB;

    int width, height;

    Spectrum[] sensor = new Spectrum[3];
    Spectrum[] ill = new Spectrum[3];
    Spectrum[] stain = new Spectrum[3];

    String[] strStain = new String[3];
    String[] strIll = new String[3];
    String[] strSensor = new String[3];

    String strMode2Stain = modes2Stain[0];
    String strModeDecon = modesDecon[0];

    boolean useNormalizedVectors = false;
    boolean usePseudo3Stain = false;
    boolean useTripleStainSeparation = false;
    boolean performDecon = true;

    String stainSubPath = "stain/";
    String illuminationSubPath = "illumination/";
    String sensorSubPath = "sensor/";

    int spectrumLength = 0;

    String spectrumFolder = "";


    public void run(String arg) {
        double[] I0;
        double[] cv1, cv2, cv3;
        double[] cv1norm=null, cv2norm=null, cv3norm=null;
        double cv1Length=1.0, cv2Length=1.0, cv3Length=1.0;

        ImagePlus imp = WindowManager.getCurrentImage();
        if (imp==null){
            IJ.error("No image!");
            return;
        }
        imp.killRoi();

        // Check processor type
        boolean grayInput = true;
        if (imp.getProcessor() instanceof ColorProcessor ){
            grayInput = false;
            impRGB = imp;
            if (impRGB.getStackSize() != 1){
                IJ.showMessage("Color image stacks are not supported");
                return;
            }
        }
        else{
            if ( !(imp.getProcessor() instanceof FloatProcessor ) ){
                if (imp.getStackSize() == 1)
                    imp.setProcessor( imp.getProcessor().convertToFloatProcessor() );
                else
                    new StackConverter(imp).convertToGray32();
            }
            if (imp.getStack().getSize() < 2 ){
                performDecon = false;
            }
        }
        width = imp.getWidth();
        height = imp.getHeight();

        // Select folder and parameters
        /*spectrumFolder = Prefs.get("dip.spetrumfolder", IJ.getDirectory("temp"));
        DirectoryChooser.setDefaultDirectory(spectrumFolder);
        DirectoryChooser dc = new DirectoryChooser("Select spectra folder");
        spectrumFolder = dc.getDirectory();*/
        spectrumFolder = "/Users/ERexhepa/Documents/TrainingTissueFinder/CDReview_Installation/CDReview_spectra/";

        if (spectrumFolder == null)
            return;
        Prefs.set("dip.spetrumfolder", spectrumFolder);

        readPrefs();
        if ( !selectParameters() )
            return;
        writePrefs();

        // Load spectra and check ..
        for (int i=0; i<3; i++){
            stain[i] = Spectrum.read(spectrumFolder + stainSubPath + strStain[i]);
            ill[i] = Spectrum.read(spectrumFolder + illuminationSubPath + strIll[i]);
            sensor[i] = Spectrum.read(spectrumFolder + sensorSubPath + strSensor[i]);
        }
        for (int i=0; i<3; i++){
            if ( ill[i] == null || sensor[i] == null ){
                IJ.showMessage("Not all necessary spectra selceted. Processing canceled!");
                return;
            }
        }

        spectrumLength = getSizeofSpectra();
        if ( spectrumLength == 0 ){
            IJ.showMessage("Spectra have different size. Check spectra files! Processing canceled");
            return;
        }

        // Calculate maximum intensity I0 (c'=0)
        I0 = new double[3];
        for (int i=0; i<spectrumLength; i++){
            I0[R] += ill[R].val[i] * sensor[R].val[i];
            I0[G] += ill[G].val[i] * sensor[G].val[i];
            I0[B] += ill[B].val[i] * sensor[B].val[i];
        }
        //IJ.log ("I0 =  " + I0[R] + "  "+ I0[G] + "  "+ I0[B]);

        // Calculate stain vectors
        cv1 = calcCV(stain[0], I0, "cv1");
        cv2 = calcCV(stain[1], I0, "cv2");
        cv3 = calcCV(stain[2], I0, "cv3");

        if ( cv1 == null ){
            IJ.showMessage("No valid stain and parameter selection. Processing canceled!");
            return;
        }
        if ( cv2 == null && cv3 == null )
            performDecon = false;
        else if (cv2 != null && cv3 != null)
            useTripleStainSeparation = true;
        else{
            if ( cv2 == null ){
                cv2 = cv3;
                stain[1] = stain[2];
                cv3 = null;
                stain[2] = null;
            }
            if ( usePseudo3Stain ){
                cv3 = calcPerpendVector(cv1.clone(), cv2.clone());
                useTripleStainSeparation = true;
            }
        }

        if ( useNormalizedVectors ){
            cv1norm = new double[3]; cv2norm = new double[3]; cv3norm = new double[3];

            cv1Length = normalizeCV(cv1, cv1norm, "cv1norm");
            IJ.log("cv1Length : " + cv1Length);
            cv2Length = normalizeCV(cv2, cv2norm, "cv2norm");
            IJ.log("cv2Length : " + cv2Length);
            cv3Length = normalizeCV(cv3, cv3norm, "cv3norm");
            IJ.log("cv3Length : " + cv3Length);

            cv1 = cv1norm; cv2 = cv2norm; cv3 = cv3norm;
        }


        // *********  Main calculations *********
        String procMode = "";
        if ( grayInput ){ // create 32bit absorbance arrays and benchmark image
            createAbsorbanceArraysAndRGBImage(imp, stain[0], stain[1], stain[2], I0);
            procMode = "_per32bit_";
            IJ.log("ready createAbsorbanceArrays");
        }
        else{ // prepare for direct deconvolution of RGB images
            prepareDirectRGBDecon(I0);
            procMode = "_per8bitRGB_";
            IJ.log("ready prepareDirectRGBDecon");
        }

        if ( performDecon ){
            if ( !useTripleStainSeparation ){
                if (strModeDecon.equals( modesDecon[0] )) // = "DET"
                    deconvolveImage2SDET(I0, cv1, cv2, procMode + modesDecon[0] + "_" + strMode2Stain);
                else if (strModeDecon.equals( modesDecon[1] )) // = "LU"
                    deconvolveImage2SLU(I0, cv1, cv2, procMode + modesDecon[1] + "_" + strMode2Stain);
                else if (strModeDecon.equals( modesDecon[2] )) // = "QR"
                    deconvolveImage2SQR(I0, cv1, cv2, procMode + modesDecon[2]);
            }
            else{
                if (strModeDecon.equals( modesDecon[0] )) // = "DET"
                    deconvolveImage3SDET(I0, cv1, cv2, cv3, procMode + modesDecon[0] );
                else if (strModeDecon.equals( modesDecon[1] )) // = "LU"
                    deconvolveImage3SLU(I0, cv1, cv2, cv3, procMode + modesDecon[1]);
                else if (strModeDecon.equals( modesDecon[2] )) // = "QR"
                    deconvolveImage3SQR(I0, cv1, cv2, cv3, procMode + modesDecon[2]);
            }
            IJ.log("ready deconvolveImage");
        }
    }


    void deconvolveImage2SDET(double[] I0, double[] cv1, double[] cv2, String ext){
        // Create Concentration Images
        impDecon = NewImage.createFloatImage("DeconResult"+ext, width, height, 1, NewImage.FILL_BLACK);
        ImageProcessor ip1 = impDecon.getProcessor();
        ImageProcessor ip2 = ip1.duplicate();
        float[] pixels1 = (float[]) ip1.getPixels();
        float[] pixels2 = (float[]) ip2.getPixels();

        int[] chn = getChannel(strMode2Stain);
        double[][] a = new double[2][];
        assignAbsorbance(a);

        double detA = cv1[chn[0]]*cv2[chn[1]] - cv1[chn[1]]*cv2[chn[0]];
        for (int i=0;i<(width*height);i++){
            pixels1[i] = (float) ((a[0][i] * cv2[chn[1]] - a[1][i] * cv2[chn[0]]) / detA);
            pixels2[i] = (float) ((cv1[chn[0]] * a[1][i] - cv1[chn[1]] * a[0][i]) / detA);
        }
        ImageStack stack = impDecon.getStack();
        stack.addSlice(ip2);
        impDecon.setStack(stack);
        impDecon.show();
    }


    void deconvolveImage3SDET(double[] I0, double[] cv1, double[] cv2, double[] cv3, String ext){
        // Create Concentration Images
        impDecon = NewImage.createFloatImage("DeconResult"+ext, width, height, 1, NewImage.FILL_BLACK);
        ImageProcessor ip1 = impDecon.getProcessor();
        ImageProcessor ip2 = ip1.duplicate();
        ImageProcessor ip3 = ip1.duplicate();
        float[] pixels1 = (float[]) ip1.getPixels();
        float[] pixels2 = (float[]) ip2.getPixels();
        float[] pixels3 = (float[]) ip3.getPixels();

        double[] a = new double[3];
        double detA = det3x3(cv1, cv2, cv3);
        for (int i=0;i<(width*height);i++){
            a[0] = aR[i]; a[1] = aG[i]; a[2] = aB[i];
            pixels1[i] = (float) ( det3x3(a, cv2, cv3) / detA);
            pixels2[i] = (float) ( det3x3(cv1, a, cv3) / detA);
            pixels3[i] = (float) ( det3x3(cv1, cv2, a) / detA);
        }
        ImageStack stack = impDecon.getStack();
        stack.addSlice(ip2);
        stack.addSlice(ip3);
        impDecon.setStack(stack);
        impDecon.show();
    }

    double det3x3(double[] v0, double[] v1, double[] v2){
        double val = v0[0]*v1[1]*v2[2] + v0[1]*v1[2]*v2[0] + v0[2]*v1[0]*v2[1]
                - v0[2]*v1[1]*v2[0] - v0[0]*v1[2]*v2[1] - v0[1]*v1[0]*v2[2];
        return val;
    }

    void deconvolveImage2SLU(double[] I0, double[] cv1, double[] cv2, String ext){
        int[] chn = getChannel(strMode2Stain);
        double[][] a = new double[2][];
        assignAbsorbance(a);

        // Prepare Linear Solver
        RealMatrix coefficients = new Array2DRowRealMatrix(new double[][] {{cv1[chn[0]], cv2[chn[0]]}, {cv1[chn[1]], cv2[chn[1]]}}, false);
        DecompositionSolver solver = new LUDecomposition(coefficients).getSolver();
        RealVector constants;
        RealVector solution;

        // Create Concentration Images
        impDecon = NewImage.createFloatImage("DeconResult"+ext, width, height, 1, NewImage.FILL_BLACK);
        ImageProcessor ip1 = impDecon.getProcessor();
        ImageProcessor ip2 = ip1.duplicate();
        float[] pixels1 = (float[]) ip1.getPixels();
        float[] pixels2 = (float[]) ip2.getPixels();

        for (int i=0;i<(width*height);i++){
            constants = new ArrayRealVector(new double[] { a[0][i], a[1][i]}, false);
            solution = solver.solve(constants);
            pixels1[i] = (float) solution.getEntry(0);
            pixels2[i] = (float) solution.getEntry(1);
        }
        ImageStack stack = impDecon.getStack();
        stack.addSlice(ip2);
        impDecon.setStack(stack);
        impDecon.show();
    }

    void deconvolveImage3SLU(double[] I0, double[] cv1, double[] cv2, double[] cv3, String ext){
        // Prepare Linear Solver
        RealMatrix coefficients = new Array2DRowRealMatrix(new double[][] {{cv1[R], cv2[R], cv3[R]}, {cv1[G], cv2[G], cv3[G]}, {cv1[B], cv2[B], cv3[B]}}, false);
        DecompositionSolver solver = new LUDecomposition(coefficients).getSolver();
        RealVector constants;
        RealVector solution;

        // Create Concentration Images
        impDecon = NewImage.createFloatImage("DeconResult"+ext, width, height, 1, NewImage.FILL_BLACK);
        ImageProcessor ip1 = impDecon.getProcessor();
        ImageProcessor ip2 = ip1.duplicate();
        ImageProcessor ip3 = ip1.duplicate();
        float[] pixels1 = (float[]) ip1.getPixels();
        float[] pixels2 = (float[]) ip2.getPixels();
        float[] pixels3 = (float[]) ip3.getPixels();

        for (int i=0;i<(width*height);i++){
            constants = new ArrayRealVector(new double[] { aR[i], aG[i], aB[i]}, false);
            solution = solver.solve(constants);
            pixels1[i] = (float) solution.getEntry(0);
            pixels2[i] = (float) solution.getEntry(1);
            pixels3[i] = (float) solution.getEntry(2);
        }
        ImageStack stack = impDecon.getStack();
        stack.addSlice(ip2);
        stack.addSlice(ip3);
        impDecon.setStack(stack);
        impDecon.show();
    }


    void deconvolveImage2SQR(double[] I0, double[] cv1, double[] cv2, String ext){
        // Prepare Linear Solver
        RealMatrix coefficients = new Array2DRowRealMatrix(new double[][] { {cv1[R], cv2[R]} , {cv1[G], cv2[G]}, {cv1[B], cv2[B]}}, false);
        DecompositionSolver solver = new QRDecomposition(coefficients).getSolver();
        RealVector constants;
        RealVector solution;

        // Create Concentration Images
        impDecon = NewImage.createFloatImage("DeconResult"+ext, width, height, 1, NewImage.FILL_BLACK);
        ImageProcessor ip1 = impDecon.getProcessor();
        ImageProcessor ip2 = ip1.duplicate();
        float[] pixels1 = (float[]) ip1.getPixels();
        float[] pixels2 = (float[]) ip2.getPixels();

        for (int i=0;i<(width*height);i++){
            constants = new ArrayRealVector(new double[] { aR[i], aG[i], aB[i]}, false);
            solution = solver.solve(constants);
            pixels1[i] = (float) solution.getEntry(0);
            pixels2[i] = (float) solution.getEntry(1);
        }
        ImageStack stack = impDecon.getStack();
        stack.addSlice(ip2);
        impDecon.setStack(stack);
        impDecon.show();
    }

    void deconvolveImage3SQR(double[] I0, double[] cv1, double[] cv2, double[] cv3, String ext){
        // Prepare Linear Solver
        RealMatrix coefficients = new Array2DRowRealMatrix(new double[][] { {cv1[R], cv2[R], cv3[R]} , {cv1[G], cv2[G], cv3[G]}, {cv1[B], cv2[B], cv3[B]}}, false);
        DecompositionSolver solver = new QRDecomposition(coefficients).getSolver();
        RealVector constants;
        RealVector solution;

        // Create Concentration Images
        impDecon = NewImage.createFloatImage("DeconResult"+ext, width, height, 1, NewImage.FILL_BLACK);
        ImageProcessor ip1 = impDecon.getProcessor();
        ImageProcessor ip2 = ip1.duplicate();
        ImageProcessor ip3 = ip1.duplicate();
        float[] pixels1 = (float[]) ip1.getPixels();
        float[] pixels2 = (float[]) ip2.getPixels();
        float[] pixels3 = (float[]) ip3.getPixels();

        for (int i=0;i<(width*height);i++){
            constants = new ArrayRealVector(new double[] { aR[i], aG[i], aB[i]}, false);
            solution = solver.solve(constants);
            pixels1[i] = (float) solution.getEntry(0);
            pixels2[i] = (float) solution.getEntry(1);
            pixels3[i] = (float) solution.getEntry(2);
        }
        ImageStack stack = impDecon.getStack();
        stack.addSlice(ip2);
        stack.addSlice(ip3);
        impDecon.setStack(stack);
        impDecon.show();
    }


    void createAbsorbanceArraysAndRGBImage(ImagePlus impC, Spectrum stain1, Spectrum stain2, Spectrum stain3, double[] I0){
        int pixelRGB;
        double transmission;
        double c1, c2, c3=0;
        int red, green, blue;
        double[] intI = new double[3];

        ImageStack stack = impC.getStack();
        int stackSize = stack.getSize();
        boolean use2Stain=false, use3Stain=false;

        float[] pixels1 = (float[]) stack.getProcessor(1).getPixels();
        float[] pixels2 = new float[pixels1.length];
        float[] pixels3 = new float[pixels1.length];
        if ( stain2 != null && stackSize >= 2 ){
            pixels2 = (float[]) stack.getProcessor(2).getPixels();
            use2Stain = true;
        }
        if ( stain3 != null && stackSize >= 3 ){
            pixels3 = (float[]) stack.getProcessor(3).getPixels();
            use3Stain = true;
        }

        // Calculate white balance values
        double maxIVal = 255.0;
        double[] wb = new double[3];
        wb[R] = maxIVal / I0[R];
        wb[G] = maxIVal / I0[G];
        wb[B] = maxIVal / I0[B];

        impRGB = IJ.createImage("RGB in", "RGB", width, height, 1);
        ImageProcessor ipRGB = impRGB.getProcessor();

        aR = new double[width*height];
        aG = new double[width*height];
        aB = new double[width*height];

        int index;
        for (int y=0; y<height; y++){
            for (int x=0; x<width; x++){
                index = y*width + x;

                c1 = pixels1[index];
                c2 = pixels2[index];
                c3 = pixels3[index];

                // Calculate transmitted/detected intensity
                intI[R] = intI[G] = intI[B] = 0.0;
                for (int i=0; i<spectrumLength; i++){
                    if ( use3Stain )
                        transmission = Math.exp(-stain1.val[i]*c1) * Math.exp(-stain2.val[i]*c2) * Math.exp(-stain3.val[i]*c3);
                    else if ( use2Stain )
                        transmission = Math.exp(-stain1.val[i]*c1) * Math.exp(-stain2.val[i]*c2);
                    else
                        transmission = Math.exp(-stain1.val[i]*c1);
                    intI[R] += ill[R].val[i] * sensor[R].val[i] * transmission;
                    intI[G] += ill[G].val[i] * sensor[G].val[i] * transmission;
                    intI[B] += ill[B].val[i] * sensor[B].val[i] * transmission;
                }

                // Calculate absorbance values
                aR[index] = -Math.log( intI[R] / I0[R] );
                aG[index] = -Math.log( intI[G] / I0[G] );
                aB[index] = -Math.log( intI[B] / I0[B] );

                // Calculate RGB image (white balanced)
                red = (int) Math.min(Math.round(intI[R] * wb[R]), 255.0);
                green = (int) Math.min(Math.round(intI[G] * wb[G]), 255.0);
                blue = (int) Math.min(Math.round(intI[B] * wb[B]), 255.0);

                pixelRGB = 0xff000000 | ((((byte)red)<<16)&0xff0000) | ((((byte)green)<<8)&0xff00) | ((byte)blue&0xff);
                ipRGB.set(x, y, pixelRGB);
            }
            IJ.showProgress(y, height);
        }
        impRGB.show();
    }


    void prepareDirectRGBDecon(double[] I0){
        // Separate RGB channels
        int[] pixels = (int[]) ((ColorProcessor)(impRGB.getProcessor())).getPixels();
        double[] pR = new double[pixels.length];
        double[] pG = new double[pixels.length];
        double[] pB = new double[pixels.length];
        int pVal;
        for (int i=0; i<pixels.length; i++){
            pVal = pixels[i];
            pR[i] = (pVal&0x00ff0000)>>16;
            pG[i] = (pVal&0x0000ff00)>>8;
            pB[i] = pVal&0x000000ff;
        }
        // Calculate absorbance values
        aR = new double[width*height];
        aG = new double[width*height];
        aB = new double[width*height];

        int index;
        for (int y=0; y<height; y++){
            for (int x=0; x<width; x++){
                index = y*width + x;
                aR[index] = -Math.log( pR[index] / 255.0 );
                aG[index] = -Math.log( pG[index] / 255.0 );
                aB[index] = -Math.log( pB[index] / 255.0 );
            }
            IJ.showProgress(y, height);
        }
    }


    double[] calcCV(Spectrum stain, double[] I0, String cvName){
        // Stain / Color Vector (c'=1)
        double [] cv;
        double [] intI = new double[3];
        cv = new double[3];

        if ( stain == null)
            return null;

        for (int i=0; i<spectrumLength; i++){
            intI[R] += ill[R].val[i] * sensor[R].val[i] * Math.exp(-stain.val[i]);
            intI[G] += ill[G].val[i] * sensor[G].val[i] * Math.exp(-stain.val[i]);
            intI[B] += ill[B].val[i] * sensor[B].val[i] * Math.exp(-stain.val[i]);
        }

        cv[R] = -(Math.log( intI[R] / I0[R] ));
        cv[G] = -(Math.log( intI[G] / I0[G] ));
        cv[B] = -(Math.log( intI[B] / I0[B] ));

        if ( !cvName.equals("") )
            IJ.log(cvName + "    " + cv[R] + "    " + cv[G] + "    " + cv[B]);

        return cv;
    }

    double normalizeCV(double[] cv, double[] cvOut, String cvName){
        if ( cv == null )
            return 0;

        double length = Math.sqrt(cv[R]*cv[R] + cv[G]*cv[G] + cv[B]*cv[B]);
        for (int i=0; i<3; i++)
            cvOut[i] = cv[i] / length;

        if ( !cvName.equals("") )
            IJ.log(cvName + "    " + cvOut[R] + "    " + cvOut[G] + "    " + cvOut[B]);

        return length;
    }


    double[] calcPerpendVector(double[] cv1, double[] cv2){
        if ( cv1 == null || cv2 == null)
            return null;

        normalizeCV(cv1, cv1, "");
        normalizeCV(cv2, cv2, "");

        double[] cvOut = new double[3];
        cvOut[R] = cv1[G]*cv2[B] - cv1[B]*cv2[G];
        cvOut[G] = cv1[B]*cv2[R] - cv1[R]*cv2[B];
        cvOut[B] = cv1[R]*cv2[G] - cv1[G]*cv2[R];

        normalizeCV(cvOut, cvOut, "");
        return cvOut;
    }

    int[] getChannel(String mode){
        int[] chn = new int[2];
        String chn1 = mode.substring(0, 1);
        String chn2 = mode.substring(1, 2);
        if ( chn1.equals("R")) chn[0]  = R;
        if ( chn1.equals("G")) chn[0]  = G;
        if ( chn1.equals("B")) chn[0]  = B;
        if ( chn2.equals("R")) chn[1]  = R;
        if ( chn2.equals("G")) chn[1]  = G;
        if ( chn2.equals("B")) chn[1]  = B;
        return chn;
    }

    void assignAbsorbance(double[][] a){
        if ( strMode2Stain.equals("GB") ){
            a[0] = aG;
            a[1] = aB;
        }
        if ( strMode2Stain.equals("GR") ){
            a[0] = aG;
            a[1] = aR;
        }
        if ( strMode2Stain.equals("RB") ){
            a[0] = aR;
            a[1] = aB;
        }
    }


    boolean selectParameters() {
        // TODO Rewrite to simplifify the dialog for simple version and totally remove for the command line version

        String[] stainSpectra = getChoices(spectrumFolder + stainSubPath);
        String[] illSpectra = getChoices(spectrumFolder + illuminationSubPath);
        String[] sensorSpectra = getChoices(spectrumFolder + sensorSubPath);

        // Reconfigure default staining parameters
        String[] stainSpectraGenfit = new String[stainSpectra.length];
        stainSpectraGenfit[0]       = "HeamatoxylinEosin";
        stainSpectraGenfit[1]       = "HeamatoxylinEosinSerafin";
        stainSpectraGenfit[2]       = "HeamatoxylinPicosirius";
        stainSpectraGenfit[3]       = "Picosirius";
        stainSpectraGenfit[4]       = "MasonTrichrome";
        stainSpectraGenfit[5]       = "HeamatoxylinDAB";

        GenericDialog gd = new GenericDialog("NASH ROBUST COLOUR DECONV - STAIN COMBINATIONS");
        gd.addChoice("STAIN PROTOCOL - ", stainSpectraGenfit, strStain[0]);
        gd.showDialog();

/*        gd.addChoice("Stain 1", stainSpectra, strStain[0]);
        gd.addChoice("Stain 2", stainSpectra, strStain[1]);
        gd.addChoice("Stain 3", stainSpectra, strStain[2]);
        gd.addChoice("Illumination R", illSpectra, strIll[R]);
        gd.addChoice("Illumination G", illSpectra, strIll[G]);
        gd.addChoice("Illumination B", illSpectra, strIll[B]);
        gd.addChoice("Sensor R", sensorSpectra, strSensor[R]);
        gd.addChoice("Sensor G", sensorSpectra, strSensor[G]);
        gd.addChoice("Sensor B", sensorSpectra, strSensor[B]);
        gd.addCheckbox("Use normalized vectors", useNormalizedVectors);
        gd.addCheckbox("Use pseudo3Stain", usePseudo3Stain);
        gd.addChoice("2 stain mode", modes2Stain, modes2Stain[0]);
        gd.addChoice("Decon mode", modesDecon, modesDecon[0]);
        gd.showDialog();*/

        if (gd.wasCanceled())
            return false;

/*        strStain[0] = gd.getNextChoice();
        strStain[1] = gd.getNextChoice();
        strStain[2] = gd.getNextChoice();
        strIll[R] = gd.getNextChoice();
        strIll[G] = gd.getNextChoice();
        strIll[B] = gd.getNextChoice();
        strSensor[R] = gd.getNextChoice();
        strSensor[G] = gd.getNextChoice();
        strSensor[B] = gd.getNextChoice();

        useNormalizedVectors = gd.getNextBoolean();
        usePseudo3Stain = gd.getNextBoolean();

        strMode2Stain = gd.getNextChoice();
        strModeDecon = gd.getNextChoice();*/
        configureStainingProtolParam(gd.getNextChoiceIndex(),stainSpectra,illSpectra,sensorSpectra);
        return true;
    }

    void configureStainingProtolParam(int choiceProtocol, String[] stainSpectra, String[] illSpectra,String[] sensorSpectra){
        switch (choiceProtocol){
            case 1:{
                strStain[0]     = stainSpectra[4];  strStain[1] = stainSpectra[2];      strStain[2] = stainSpectra[0];
                strIll[R]       = illSpectra[2];    strIll[G]   = illSpectra[2];        strIll[B]   = illSpectra[2];
                strSensor[R]    = strSensor[4];     strSensor[G]= strSensor[3];         strSensor[B]= strSensor[1];

                useNormalizedVectors    = Boolean.TRUE;
                usePseudo3Stain         = Boolean.FALSE;

                strMode2Stain           = this.modes2Stain[2];
                strModeDecon            = this.modesDecon[2];
            }
            case 2:{
                strStain[0]     = stainSpectra[4];  strStain[1] = stainSpectra[2];      strStain[2] = stainSpectra[0];
                strIll[R]       = illSpectra[2];    strIll[G]   = illSpectra[2];        strIll[B]   = illSpectra[2];
                strSensor[R]    = sensorSpectra[4];     strSensor[G]= sensorSpectra[3];         strSensor[B]= sensorSpectra[1];

                useNormalizedVectors    = Boolean.TRUE;
                usePseudo3Stain         = Boolean.FALSE;

                strMode2Stain           = this.modes2Stain[2];
                strModeDecon            = this.modesDecon[2];
            }
            case 3:{
                strStain[0]     = stainSpectra[4];  strStain[1] = stainSpectra[2];      strStain[2] = stainSpectra[0];
                strIll[R]       = illSpectra[2];    strIll[G]   = illSpectra[2];        strIll[B]   = illSpectra[2];
                strSensor[R]    = sensorSpectra[4];     strSensor[G]= sensorSpectra[3];         strSensor[B]= sensorSpectra[1];

                useNormalizedVectors    = Boolean.TRUE;
                usePseudo3Stain         = Boolean.FALSE;

                strMode2Stain           = this.modes2Stain[2];
                strModeDecon            = this.modesDecon[2];
            }
            case 4:{
                strStain[0]     = stainSpectra[4];  strStain[1] = stainSpectra[2];      strStain[2] = stainSpectra[0];
                strIll[R]       = illSpectra[2];    strIll[G]   = illSpectra[2];        strIll[B]   = illSpectra[2];
                strSensor[R]    = sensorSpectra[4];     strSensor[G]= sensorSpectra[3];         strSensor[B]= sensorSpectra[1];

                useNormalizedVectors    = Boolean.TRUE;
                usePseudo3Stain         = Boolean.FALSE;

                strMode2Stain           = this.modes2Stain[2];
                strModeDecon            = this.modesDecon[2];
            }
            case 5:{
                strStain[0]     = stainSpectra[4];  strStain[1] = stainSpectra[2];      strStain[2] = stainSpectra[0];
                strIll[R]       = illSpectra[2];    strIll[G]   = illSpectra[2];        strIll[B]   = illSpectra[2];
                strSensor[R]    = sensorSpectra[4];     strSensor[G]= sensorSpectra[3];         strSensor[B]= sensorSpectra[1];

                useNormalizedVectors    = Boolean.TRUE;
                usePseudo3Stain         = Boolean.FALSE;

                strMode2Stain           = this.modes2Stain[2];
                strModeDecon            = this.modesDecon[2];
            }
            default:{
                // Default case is HE
                strStain[0]     = stainSpectra[4];  strStain[1] = stainSpectra[2];      strStain[2] = stainSpectra[0];
                strIll[R]       = illSpectra[2];    strIll[G]   = illSpectra[2];        strIll[B]   = illSpectra[2];
                strSensor[R]    = sensorSpectra[4];     strSensor[G]= sensorSpectra[3];         strSensor[B]= sensorSpectra[1];

                useNormalizedVectors    = Boolean.TRUE;
                usePseudo3Stain         = Boolean.FALSE;

                strMode2Stain           = this.modes2Stain[2];
                strModeDecon            = this.modesDecon[2];
            }
        }

    }

    void readPrefs(){
        strStain[0] = Prefs.get("dip.strStain1", "");
        strStain[1] = Prefs.get("dip.strStain2", "");
        strStain[2] = Prefs.get("dip.strStain3", "");
        strIll[R] = Prefs.get("dip.strIllChnR", "");
        strIll[G] = Prefs.get("dip.strIllChnG", "");
        strIll[B] = Prefs.get("dip.strIllChnB", "");
        strSensor[R] = Prefs.get("dip.strSensorChnR", "");
        strSensor[G] = Prefs.get("dip.strSensorChnG", "");
        strSensor[B] = Prefs.get("dip.strSensorChnB", "");
    }

    void writePrefs(){
        Prefs.set("dip.strStain1", strStain[0]);
        Prefs.set("dip.strStain2", strStain[1]);
        Prefs.set("dip.strStain3", strStain[2]);
        Prefs.set("dip.strIllChnR", strIll[R]);
        Prefs.set("dip.strIllChnG", strIll[G]);
        Prefs.set("dip.strIllChnB", strIll[B]);
        Prefs.set("dip.strSensorChnR", strSensor[R]);
        Prefs.set("dip.strSensorChnG", strSensor[G]);
        Prefs.set("dip.strSensorChnB", strSensor[B]);
    }


    String[] getChoices(String path){
        ArrayList<String> list = new ArrayList<String>();
        String name;
        File fPath = new File(path);
        File[] files = fPath.listFiles();
        for (int i=0; i<files.length; i++){
            name = files[i].getName();
            if ( files[i].isFile() && name.endsWith(".txt") )
                list.add(name);
        }
        String[] choices = new String[list.size()+1];
        choices[0] = "--";
        for (int i=0; i<list.size(); i++)
            choices[i+1] = list.get(i);

        return choices;
    }

    int getSizeofSpectra(){
        int specSize = ill[R].getSize();

        if (specSize < 1)
            return 0;

        if ( stain[0] != null && stain[0].getSize() != specSize) return 0;
        if ( stain[1] != null && stain[1].getSize() != specSize) return 0;
        if ( stain[2] != null && stain[2].getSize() != specSize) return 0;

        if ( ill[R] != null && ill[R].getSize() != specSize) return 0;
        if ( ill[G] != null && ill[G].getSize() != specSize) return 0;
        if ( ill[B] != null && ill[B].getSize() != specSize) return 0;

        if ( sensor[R] != null && sensor[R].getSize() != specSize) return 0;
        if ( sensor[G] != null && sensor[G].getSize() != specSize) return 0;
        if ( sensor[B] != null && sensor[B].getSize() != specSize) return 0;

        return specSize;
    }

}

