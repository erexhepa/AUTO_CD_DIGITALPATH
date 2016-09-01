package AUTO_COL_DECONV;

/***************************************************
 * 
 * License
 * 
 * PlugIn:		CDReview
 * 
 * File: 		CDR_CreatePatternImage.java
 * 
 * Description:	CDReview CDR_CreatePatternImage - Software for 
 * 				creating 32bit channel patter image stack
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




import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.NewImage;
import ij.plugin.LutLoader;
import ij.plugin.PlugIn;
import ij.process.ImageProcessor;

public class CDR_CreatePatternImage implements PlugIn {

public void run(String arg) {

	int nElem = 10;
	float[] val1 = {0.0f, 0.125f, 0.25f, 0.5f, 1.0f, 1.5f, 2.0f, 3.0f, 4.0f, 5.0f};
	float[] val2 = {0.0f, 0.125f, 0.25f, 0.5f, 1.0f, 1.5f, 2.0f, 3.0f, 4.0f, 5.0f};
	float[] val3 = {0.0f, 0.125f, 0.25f, 0.5f, 1.0f, 1.5f, 2.0f, 3.0f, 4.0f, 5.0f};
	
	float val1mult = 4.0f/ 4.0f;
	float val2mult = 4.0f/ 4.0f;
	float val3mult = 4.0f/ 4.0f;
	
	int r3Width = 10;
	int r3Spacing= 5;
	int r2Spacing = 10;
	int r1Spacing = 30;
	
	int r3Height = r3Width;
	int r2Width = r3Width + 2 * r3Spacing;
	int r2Height = nElem * r3Height + (nElem+1) * r3Spacing;
	int r1Width = nElem * r2Width + (nElem+1) * r2Spacing;
	int r1Height = r2Height + + 2 * r2Spacing;
	
	int width = nElem * r1Width + (nElem+1) * r1Spacing;
	int height = r1Height + + 2 * r1Spacing;
	
	ImagePlus imp = NewImage.createFloatImage("3Stain_Concentration", width, height, 1, NewImage.FILL_BLACK);
	ImageProcessor ip1 = imp.getProcessor();
	ImageProcessor ip2 = ip1.duplicate();
	ImageProcessor ip3 = ip1.duplicate();
	
	for (int i=0; i<nElem; i++){
		ip1.setRoi((i+1)*r1Spacing + i*r1Width, r1Spacing, r1Width, r1Height);
		ip1.setColor(val1[i] * val1mult);
		ip1.fill();		
	}
	
	for (int i=0; i<nElem; i++){
		for (int k=0; k<nElem; k++){
			ip2.setRoi((i+1)*r1Spacing + i*r1Width + (k+1)*r2Spacing + k*r2Width, r1Spacing+r2Spacing, r2Width, r2Height);
			ip2.setColor(val2[k] * val2mult);
			ip2.fill();
		}
	}
	
	for (int i=0; i<nElem; i++){
		for (int k=0; k<nElem; k++){
			for (int l=0; l<nElem; l++){
				ip3.setRoi((i+1)*r1Spacing + i*r1Width + (k+1)*r2Spacing + k*r2Width + r3Spacing, r1Spacing+r2Spacing + (l+1)*r3Spacing + l*r3Height, r3Width, r3Height);
				ip3.setColor(val3[l] * val3mult);
				ip3.fill();
			}
		}
	}
		
	ip1.resetMinAndMax();
	ImageStack stack = imp.getStack();
	stack.addSlice(ip2);
	stack.addSlice(ip3);
	imp.setStack(stack);
 	imp.show();
    
	LutLoader lutLoader = new LutLoader();
	lutLoader.run("fire");

}

}
