//scaling factors and variables: th and s-th for example
standarddev = 1.5;
shaftstandarddev = 0.65;
//Imaging analysis program
run("Close All");
if (isOpen("ROI Manager")) {
     selectWindow("ROI Manager");
     run("Close");
  }
   if (isOpen("Results")) {
     selectWindow("Results");
     run("Close");
  }
  
Dialog.create("Timepoints");
  Dialog.addNumber("n:", 2);
  Dialog.show();
  n = Dialog.getNumber();
  m= n+1;
open();
path = getDirectory("image");
name = File.nameWithoutExtension;
run("Z Project...", "projection=[Max Intensity] all");
rename("image");

waitForUser("Select a region of dendrite then crop. press OK");
run("Subtract Background...", "rolling=50 stack");
// Split the colors, smooth and rename them
run("Split Channels");

selectWindow("C4-image");
	rename("fill");
selectWindow("C2-image");
	rename("camk");
selectWindow("C3-image");
	runMacro("EnhancePuncta");
	rename("psd");
selectWindow("C1-image");
runMacro("EnhancePuncta");
	rename("geph");
run("Set Measurements...", "mean standard limit redirect=None decimal=2");

//threshold whole cell with whichever channel is best
selectWindow("fill");
		run("Bleach Correction", "correction=[Histogram Matching]");
		run("Gaussian Blur...", "sigma=1 stack");
		rename("cell");
		/*run("Threshold...");
			waitForUser("select entire cell, press OK");
		run("Close");*/
		setAutoThreshold("Otsu dark stack");
		run("Convert to Mask", "method=Otsu background=Dark calculate");
	selectWindow("cell");
		run("Duplicate...", "duplicate");
			rename("t_cell");
				run("Divide...", "value=255 stack");
				setMinAndMax(0, 0);

//determine psd95 threshold
imageCalculator("Multiply create stack", "psd","t_cell");
setThreshold(1, 65535);
run("Measure");
mean = getResult("Mean",0);
sd = getResult("StdDev",0);
th = mean + ((standarddev)*sd);
s_th = mean + ((shaftstandarddev)*sd);
selectWindow("Results");
run("Close");

//CreateImagesforThreshold
selectWindow("psd");
	run("Duplicate...", "duplicate");
		rename("puncta");
			run("Duplicate...", "duplicate");
				rename("p_shaft");
selectWindow("puncta");
//run("Threshold...");
	setAutoThreshold("Default dark no-reset");
		setThreshold(th, 65535);
			setOption("BlackBackground", false);
			run("Convert to Mask", "method=Default background=Dark");
			
selectWindow("p_shaft");
	setAutoThreshold("Default dark no-reset");
//run("Threshold...");
		setThreshold(s_th, 65535);
			setOption("BlackBackground", false);
			run("Convert to Mask", "method=Default background=Dark");
				run("Dilate", "stack");
				
//Determine geph threshold
imageCalculator("Multiply create stack", "geph","t_cell");
setThreshold(1, 65535);
run("Measure");
mean = getResult("Mean",0);
sd = getResult("StdDev",0);
th = mean + ((standarddev)*sd);
s_th = mean + ((shaftstandarddev)*sd);
selectWindow("Results");
run("Close");

selectWindow("geph");
	run("Duplicate...", "duplicate");
		rename("g_puncta");
			run("Duplicate...", "duplicate");
				rename("g_shaft");
selectWindow("g_puncta");
	//run("Threshold...");
	setAutoThreshold("Default dark no-reset");
		setThreshold(th, 65535);
			setOption("BlackBackground", false);
			run("Convert to Mask", "method=Default background=Dark");
	
selectWindow("g_shaft");
	setAutoThreshold("Default dark no-reset");
//run("Threshold...");
		setThreshold(s_th, 65535);
			setOption("BlackBackground", false);
			run("Convert to Mask", "method=Default background=Dark");
				run("Dilate", "stack");

//create actual shaft mask. Subtracting dilated "psd" mask (shaft1) from whole cell leaves just the shaft!
imageCalculator("Subtract create stack", "cell","p_shaft");
	rename("shaftPSD");
	imageCalculator("Subtract create stack", "shaftPSD","g_shaft");
	rename("shaft");
		run("Open", "stack");
		setOption("BlackBackground", false);
		run("Erode", "stack");
		run("Divide...", "value=255 stack");
		setMinAndMax(0, 0);

selectWindow("shaft");
run("Duplicate...", "title=ShaftStack duplicate");
//here we go! this is the analysis stage. Now we can measure actual images with thesholded objects to measure average intensity values within objects
run("Clear Results");
//Overlay the shaft threshold onto the color of interest
imageCalculator("Multiply create stack", "camk","shaft");

// Measure sum of all pixel intensities from p_shaft mask of color of interest
run("Set Measurements...", "mean limit redirect=None decimal=2");
run("Stack to Images");
for (i=1; i<m; i++) {
	selectWindow("Result-000"+i);
		setThreshold(1, 65535);
		run("Measure");
}

// measurements here are the number of pixels from each shaft mask. This is the area.

//turning shaft values into an array for computation
shaftvalues = newArray(getResult("Mean",0), getResult("Mean",1));
for (j=2; j<nResults; j++) {
	shaftvalues = Array.concat(shaftvalues, getResult("Mean",j));
}
run("Clear Results");
	
//Now is when we take the average intensity of CaMKII within each PSD object. The data that results will be 1) the area of each object and the centroid coordinates, and 2) the average CaMKII intensity within each object.
//This program will cycle through the timepoints until no timepoints remain
selectWindow("puncta");
run("Duplicate...", "duplicate");
rename("punctastack");
selectWindow("puncta");
run("Stack to Images");
selectWindow("camk");
run("Stack to Images");
selectWindow("psd");
run("Stack to Images");
selectWindow("g_puncta");
run("Duplicate...", "duplicate");
rename("gephstack");
selectWindow("g_puncta");
run("Stack to Images");

for (i=1; i<m; i++) {
	selectWindow("puncta-000"+i);
	run("Set Measurements...", "area centroid redirect=None decimal=2");
	run("Analyze Particles...", "size=2-Infinity pixel show=Outlines display add");

	PSDlength = nResults;
	objectarea = newArray(getResult("Area",0), getResult("Area",1));
	for (j=2; j<nResults; j++) {
		objectarea = Array.concat(objectarea, getResult("Area",j));
	}
	xcentroids = newArray(getResult("X",0), getResult("X",1));
	for (j=2; j<nResults; j++) {
		xcentroids = Array.concat(xcentroids, getResult("X",j));
	}
	ycentroids = newArray(getResult("Y",0), getResult("Y",1));
	for (j=2; j<nResults; j++) {
		ycentroids = Array.concat(ycentroids, getResult("Y",j));
	}
				run("Clear Results");
		run("Set Measurements...", "area mean redirect=None decimal=2");
		selectWindow("camk-000"+i);
			roiManager("multi-measure measure_all");
	//an array of camkii values that is the average intensity within the object divided by the mean intensity of camkii in the shaft (intensity/num pixels)
	
	camkiiPSD = newArray(
		getResult("Mean",0)/(shaftvalues[i-1]),
		getResult("Mean",1)/(shaftvalues[i-1]));
	for (j=2; j<PSDlength; j++) {
		camkiiPSD = Array.concat(camkiiPSD,
			getResult("Mean",j)/(shaftvalues[i-1]));
	}
		run("Clear Results");


	run("Clear Results");
	selectWindow("ROI Manager");
	run("Close");
	//geph puncta centroids
		run("Set Measurements...", "area centroid redirect=None decimal=2");
		selectWindow("g_puncta-000"+i);
		run("Analyze Particles...", "size=2-Infinity pixel show=Outlines display add");
	gephobjectarea = newArray(getResult("Area",0), getResult("Area",1));
	for (j=2; j<nResults; j++) {
		gephobjectarea = Array.concat(gephobjectarea, getResult("Area",j));
	}
	gephxcentroids = newArray(getResult("X",0), getResult("X",1));
	for (j=2; j<nResults; j++) {
		gephxcentroids = Array.concat(gephxcentroids, getResult("X",j));
	}
	gephycentroids = newArray(getResult("Y",0), getResult("Y",1));
	for (j=2; j<nResults; j++) {
		gephycentroids = Array.concat(gephycentroids, getResult("Y",j));
	}
	gephlength = nResults;
	run("Clear Results");			
		run("Set Measurements...", "area mean redirect=None decimal=2");
		selectWindow("camk-000"+i);
			roiManager("multi-measure measure_all");
	//an array of camkii values that is the average intensity within the object divided by the mean intensity of camkii in the shaft (intensity/num pixels)
	
	camkiiGeph = newArray(
		getResult("Mean",0)/(shaftvalues[i-1]),
		getResult("Mean",1)/(shaftvalues[i-1]));
	for (j=2; j<gephlength; j++) {
		camkiiGeph = Array.concat(camkiiGeph,
			getResult("Mean",j)/(shaftvalues[i-1]));
	}
	run("Clear Results");
	selectWindow("ROI Manager");
	run("Close");
	
	minz = newArray(9999,9999);
	for(j=0; j<PSDlength; j++) {
		Zei = newArray(9999,9999);
		Ex = xcentroids[j];
		Ey = ycentroids[j];
	
		for (k=0; k<gephlength; k++) {
			Ix = gephxcentroids[k];
			Iy = gephycentroids[k];
				z =(sqrt( ((Ex-Ix)*(Ex-Ix)) + ((Ey-Iy)*(Ey-Iy))));
					Zei = Array.concat(Zei, z );
			}
	rank = Array.rankPositions(Zei);
	minz = Array.concat(minz, Zei[rank[0]]);
	
	}
	minz = Array.deleteIndex(minz, 0);
	minz = Array.deleteIndex(minz, 0);	
	run("Clear Results");

//psd proximity to geph
	g_minz = newArray(9999,9999);
	for (k=0; k<gephlength; k++) {
			Zei = newArray(9999,9999);
			Ix = gephxcentroids[k];
			Iy = gephycentroids[k];
			
		for(j=0; j<PSDlength; j++) {
			Ex = xcentroids[j];
			Ey = ycentroids[j];
				z =(sqrt( ((Ex-Ix)*(Ex-Ix)) + ((Ey-Iy)*(Ey-Iy))));
					Zei = Array.concat(Zei, z );
			}
	rank = Array.rankPositions(Zei);
	g_minz = Array.concat(g_minz, Zei[rank[0]]);
	
	}
	g_minz = Array.deleteIndex(g_minz, 0);
	g_minz = Array.deleteIndex(g_minz, 0);
	Array.show("Spine/shaft ratios",minz,camkiiPSD,objectarea,g_minz,camkiiGeph,gephobjectarea);
		saveAs("Results",path+"data"+name+"_timepoint_"+i+".csv");
}
selectWindow("Results");
run("Close");
run("Merge Channels...", "c1=ShaftStack c2=punctastack c3=gephstack create keep");
//run("Brightness/Contrast...");
run("Enhance Contrast", "saturated=0.35");
saveAs("tiff",path+"Masks"+name);
waitForUser("Everything look OK?");
run("Close All");

