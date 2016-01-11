

//   /home/mmeuli/.imagej/IJ_Prefs.txt
//   ij.y=14
//   ij.x=92

//E data: e_PC
//E length: Sphere_Regularization
//Initializatioin: LocalMax
//Region Competition Parameters
//no fusion, no fission, no handles
//Lambda E length: 0.0400
//Thete E merge: 0.000
//Max Iterations: 300
//Oscillation threshold (convergence): 0.0200
//Curvature based gradient flow option: 8
//Local Max Initialization:
//Radius: 10
//Sigma: 3
//Tolerance: 0.5
//Region Tol: 20

pixelW=0.04505
pixelH=0.04505
pixelD=0.15020

x=100; y=130;
call("ij.gui.WaitForUserDialog.setNextLocation",x,y);
waitForUser("Run Region competition before this Macro to set all parameters correct");

//  inputDir = getDirectory("Choose the Directory of the input files. ");
inputDir = "/media/mmeuli/WD-HD-ext4/20150903_BCG_Pasteur-Aeras_zmp1_ko_in_RAW/20160104-ijm-h5/in/"
//  inputDir = getDirectory("Choose the Directory for the output files. ");
outputDir = "/media/mmeuli/WD-HD-ext4/20150903_BCG_Pasteur-Aeras_zmp1_ko_in_RAW/20160104-ijm-h5/out/"
list = getFileList(inputDir);
count = 0;
for (i=0; i<list.length; i++) {
	if (endsWith(list[i], ".h5")) {
		count++;
	}
}

print("Number of .h5 files = " + count);
h5files = newArray(count);
h5count = 0;
for (i=0; i<list.length; i++) {
	if (endsWith(list[i], ".h5")) {
		h5files[h5count] = inputDir + list[i];
		print((h5count) + ": " + list[i]);
		h5count++;
	}
}

script = 
"lw = WindowManager.getFrame('Log');\n"+ 
"if (lw!=null) {\n"+ 
"   lw.setLocation(100,250);\n"+ 
"   lw.setSize(800, 900)\n"+ 
"}\n"; 
eval("script", script); 
	  
call("ij.gui.WaitForUserDialog.setNextLocation",x,y);
waitForUser("Check filenumber to start with:");
Dialog.create("Enter filenumber:");
Dialog.addNumber("Filenumber:", 0);
Dialog.setLocation(x,y);
Dialog.show();
filenumber = Dialog.getNumber();

while (filenumber < h5count) {	
	filename = h5files[filenumber];
	fileT = File.getName(filename); 
	columns=split(fileT,"-"); 
	cNr=parseInt(columns[0]); 
	col1=columns[1]; 
	dotIndex = indexOf(col1, "."); 
	iNr = substring(col1, 0, dotIndex); 

	print("Opening file: " + filename);
	run("Load HDF5 File...", "open=filename");
	call("ij.gui.WaitForUserDialog.setNextLocation",x,y);
	waitForUser("Load data set of h5 file and then click ok.");
	title = getTitle(); 
	selectWindow(title);
	run("Set... ", "zoom=150"); 
	setLocation(100, 220);
	run("Properties...", "pixel_width=pixelW pixel_height=pixelH voxel_depth=pixelD");
	getDimensions(w, h, channels, slices, frames);
	m_slice = slices/2;
	Stack.setPosition(1, m_slice, 1);
	run("Green"); 
	run("Enhance Contrast", "saturated=0.35");
	Stack.setPosition(2, m_slice, 1);
	run("Red"); 
//	setMinAndMax(0, 20);
	run("Enhance Contrast", "saturated=0.35");
	Stack.setPosition(3, m_slice, 1);
	run("Blue"); 
	run("Enhance Contrast", "saturated=0.35");
	wait(500);
	for (i=1; i<=slices; i++) {
	Stack.setPosition(1, i, 1);
	wait(100);
	}
	Stack.setPosition(1, m_slice, 1);

	call("ij.gui.WaitForUserDialog.setNextLocation",x,y);
	waitForUser("If there are bacteria(s) draw a rectangular to crop the image.\nIf not deselect checkbox after this dialog.:");
	Dialog.create("Bacteria selection");
	Dialog.setLocation(x,y);
	Dialog.addCheckbox("There are bacterias on the image", true);
	Dialog.show();
	bac_check = Dialog.getCheckbox();
	if (bac_check==true) {
		run("Crop");
		run("Split Channels");
		c1Title = "C1-" + title; 
		c2Title = "C2-" + title; 
		c3Title = "C3-" + title; 
		selectWindow(c1Title);
		run("Set... ", "zoom=200"); 
		setLocation(100, 220);
		selectWindow(c2Title);
		run("Set... ", "zoom=200"); 
		setLocation(500, 220);
		selectWindow(c3Title);
		run("Set... ", "zoom=200"); 
		setLocation(100, 620);
		
		selectWindow(c1Title);
		run("Region Competition", "e_data=e_PC e_length=Sphere_Regularization lambda=0.0400 theta=0.0000 max_iterations=300 oscillation=0.0200 initialization=LocalMax inputimage=[c1Title] labelimage=[] keep_frames show_and_save_statistics");
		titleSeg = getTitle(); 
		selectWindow(titleSeg);
		run("Set... ", "zoom=200"); 
		setLocation(500, 620);
		run("Properties...", "pixel_width=pixelW pixel_height=pixelH voxel_depth=pixelD");
		Stack.setSlice(m_slice); 
		run("Enhance Contrast", "saturated=0.35"); 
	
	
		//f = File.open(""); // display file open dialog
		//Format:	cNr	iNr	label	x	y	z	pSize	pixels	maxDiameter	minDiameter	roundness	lysosome	macrophage
		resultF = outputDir + cNr + "-" + iNr + "-" + "result.txt";
		f = File.open(resultF);
	
		run("3D Manager");
		selectWindow(titleSeg);
		Ext.Manager3D_AddImage();
		call("ij.gui.WaitForUserDialog.setNextLocation",x,y);
		waitForUser("Deselect unwanted objects.");
	
		Ext.Manager3D_Count(nb_obj);
		print("number of objects",nb_obj);
	
		for (object=0; object < nb_obj; object++) {
			selectWindow(c2Title);
			//The object number and the type of measure ("IntDen", "Mean", "Min", "Max", "Sigma") 
			Ext.Manager3D_Quantif3D(object, "Mean", mean);
			print("Mean of object in C2" + object + " = " + mean);
			Ext.Manager3D_Centroid3D(object,cx,cy, cz);
			//The object number and the type of measure ("Vol", "Surf", “NbVox”, "Comp", "Feret", "Elon1", "Elon2", "DCMin", "DCMax", "DCMean", "DCSD")	
			Ext.Manager3D_Measure3D(object,"Vol",volume);	
			Ext.Manager3D_Measure3D(object,"NbVox",nbVox);
			Ext.Manager3D_Measure3D(object,"Elon1",elon1);
			Ext.Manager3D_Measure3D(object,"Elon2",elon2);
			selectWindow(c3Title);
			//The object number and the type of measure ("IntDen", "Mean", "Min", "Max", "Sigma") 
			Ext.Manager3D_Quantif3D(object, "Mean", meanM);
			print("Mean of object in C3" + object + " = " + meanM);
			print(f, cNr + "\t" + iNr + "\t" + object + "\t" + (cx*pixelW) + "\t" + (cy*pixelH) + "\t" + (cz*pixelD) + "\t" + volume + "\t" + nbVox + "\t" + elon1 + "\t" + elon2 + "\t" + 1 + "\t" + d2s(mean,3) + "\t" + meanM + "\n");
		}
		Ext.Manager3D_Close();
		File.close(f);
	
		selectWindow(c1Title);
		resultB = outputDir + cNr + "-" + iNr + "-" + "bacteria.tif";
		saveAs("Tiff", resultB);
		selectWindow(c2Title);
		resultL = outputDir + cNr + "-" + iNr + "-" + "lysosomes.tif";
		saveAs("Tiff", resultL);
		selectWindow(c3Title);
		resultM = outputDir + cNr + "-" + iNr + "-" + "macrophage.tif";
		saveAs("Tiff", resultM);
		selectWindow(titleSeg);
		resultS = outputDir + cNr + "-" + iNr + "-" + "labeled.tif";
		saveAs("Tiff", resultS);
	
		call("ij.gui.WaitForUserDialog.setNextLocation",x,y);
		waitForUser("Check results and click ok to continue.");
	} 
	else {
		print ("Image" + filename + "skiped, because there were no bacteria");
	}		



	while (nImages>0) { 
		selectImage(nImages); 
		close(); 
	}
	list = getList("window.titles"); 
	for (i=0; i<list.length; i++){ 
		winame = list[i]; 
		print(winame);
		selectWindow(winame); 
		run("Close"); 
	} 
	filenumber++;
}





 
  
