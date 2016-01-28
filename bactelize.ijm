//setBatchMode(true);
x=1600; y=900;
minVol = 0.6;	// excluded objects will in a green circle in the composite image
maxVol = 5; 		// excluded objects will  be in a white circle in the composite image 
minMacroValue = 400;   	// excluded objects will be in a red circle in the composite image
minMacroValueRed = 500;  	// objects below value will be red in the RGB macrophage image, otherwise objects are coloured from green to white according to size

//  inputDir = getDirectory("Choose the Directory of the input files. ");
inputDir = "/media/mmeuli/WD-HD-ext4/20150325_BCG_Pasteur-Aeras_zmp1_ko_in_RAW/20160128/in/"
//  inputDir = getDirectory("Choose the Directory for the output files. ");
outputDir = "/media/mmeuli/WD-HD-ext4/20150325_BCG_Pasteur-Aeras_zmp1_ko_in_RAW/20160128/out/"

list = getFileList(outputDir);
count = 0;
for (i=0; i<list.length; i++) {
	if (endsWith(list[i], "-labeled.tif")) {
		count++;
	}
}
print("Number of -labeled.tif files: " + count);
segFiles = newArray(count);
segFilesCount = 0;
for (i=0; i<list.length; i++) {
	if (endsWith(list[i], "-labeled.tif")) {
		segFiles[segFilesCount] = outputDir + list[i];
		print((segFilesCount) + ": " + list[i]);
		segFilesCount++;
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

while (filenumber < segFilesCount) {
	repeatcomposite	= false;
	filename = segFiles[filenumber];
	fileT = File.getName(filename); 
	columns=split(fileT,"-"); 
	cNr=parseInt(columns[0]); 
	iNr = parseInt(columns[1]); 
	print("Opening file: " + filename);
	open(filename);
	titleSeg = getTitle(); 
	selectWindow(titleSeg);
	run("Set... ", "zoom=150"); 
	setLocation(100, 220);
	getDimensions(w, h, channels, slices, frames);
	m_slice = slices/2;
	Stack.setPosition(1, m_slice, 1);
	run("Clear Results");
	run("Statistics"); 
	min=getResult("Min", 0);
	max=getResult("Max", 0);
	selectWindow("Results"); 
	run("Close");
	setMinAndMax(min, max);
	run("3-3-2 RGB");
	for (i=1; i<=slices; i++) {
		Stack.setPosition(1, i, 1);
		wait(50);
	}
	for (i=slices; i>=1; i--) {
		Stack.setPosition(1, i, 1);
		wait(50);
	}	
	Stack.setPosition(1, m_slice, 1);
	
	open(inputDir + cNr + "-" + iNr + "-lysosomes.tif");
	c2Title = getTitle(); 
	selectWindow(c2Title);
	run("Set... ", "zoom=70");
	setLocation(1000, 80);
	setSlice(nSlices/2);
	run("Red"); 
	run("Enhance Contrast", "saturated=1");
	
	open(inputDir + cNr + "-" + iNr + "-macrophage.tif");
	c3Title = getTitle(); 
	selectWindow(c3Title);	
	run("Set... ", "zoom=70");
	setLocation(1500, 80);
	setSlice(nSlices/2);
	run("Blue"); 
	run("Enhance Contrast", "saturated=1");

	selectWindow(c3Title);	
	getVoxelSize(width, height, depth, unit);
	selectWindow(titleSeg);
	run("Properties...", "pixel_width=width pixel_height=height voxel_depth=depth");
	
//	run("3D Manager Options", "exclude_objects_on_edges_xy exclude_objects_on_edges_z distance_between_centers=10");    
	run("3D Manager");
	selectWindow(titleSeg);
	Ext.Manager3D_AddImage();
	Ext.Manager3D_MultiSelect();
	Ext.Manager3D_Count(nb_obj);
	logF = outputDir + cNr + "-" + iNr + "-" + "excluded.txt";
	lf = File.open(logF);
	for (object=0; object < nb_obj; object++) {
		selectWindow(c2Title);
		//The object number and the type of measure ("IntDen", "Mean", "Min", "Max", "Sigma") 
		Ext.Manager3D_Quantif3D(object, "Mean", mean);;
		Ext.Manager3D_Centroid3D(object, cx, cy, cz);
		//The object number and the type of measure ("Vol", "Surf", “NbVox”, "Comp", "Feret", "Elon1", "Elon2", "DCMin", "DCMax", "DCMean", "DCSD")	
		Ext.Manager3D_Measure3D(object,"Vol",volume);
		Ext.Manager3D_Measure3D(object,"NbVox",nbVox);
		Ext.Manager3D_Measure3D(object,"Elon1",elon1);
		Ext.Manager3D_Measure3D(object,"Elon2",elon2);
		selectWindow(c3Title);
		//The object number and the type of measure ("IntDen", "Mean", "Min", "Max", "Sigma") 
		Ext.Manager3D_Quantif3D(object, "Mean", meanM);
		selectWindow(titleSeg);
		Ext.Manager3D_Quantif3D(object, "Mean", meanSeg);
		if (volume > maxVol) {  
			Ext.Manager3D_Select(object);
			drawEllipseL((cx-30), (cy-30), "white", 2);
			print(lf, cNr + "\t" + iNr + "\t" + meanSeg + "\t" + (cx*width) + "\t" + (cy*height) + "\t" + (cz*depth) + "\t" + volume + "\t" + nbVox + "\t" + elon1 + "\t" + elon2 + "\t" + 1 + "\t" + d2s(mean,3) + "\t" + meanM + "\n");
		}
		else if (volume < minVol) {  
			Ext.Manager3D_Select(object);
			drawEllipseS((cx-20), (cy-20), "green", 2);
			print(lf, cNr + "\t" + iNr + "\t" + meanSeg + "\t" + (cx*width) + "\t" + (cy*height) + "\t" + (cz*depth) + "\t" + volume + "\t" + nbVox + "\t" + elon1 + "\t" + elon2 + "\t" + 1 + "\t" + d2s(mean,3) + "\t" + meanM + "\n");
		}
		else if (meanM < minMacroValue) {  
			Ext.Manager3D_Select(object);
			drawEllipseM((cx-25), (cy-25), "red", 2);  
			print(lf, cNr + "\t" + iNr + "\t" + meanSeg + "\t" + (cx*width) + "\t" + (cy*height) + "\t" + (cz*depth) + "\t" + volume + "\t" + nbVox + "\t" + elon1 + "\t" + elon2 + "\t" + 1 + "\t" + d2s(mean,3) + "\t" + meanM + "\n");
			repeatcomposite	= true;
		}
	}
	File.close(lf);
	selectWindow(titleSeg);
	Ext.Manager3D_Delete();	

	selectWindow(c3Title);
	Ext.Manager3D_Count(nb);
	Ext.Manager3D_Measure3D(0,"Vol",V);
	Ext.Manager3D_Quantif3D(object, "Mean", meanM);
	max=V; min=V;
	// loop to find max and min volumes
	for(i=1;i<nb;i++) {
		Ext.Manager3D_Measure3D(i,"Vol",V);
		if(V>max) {
			max=V;
		}
		if(V<min){
			min=V;
		}		
	} 
	run("Blue"); 
	run("Enhance Contrast", "saturated=1");
	run("Duplicate...", "duplicate");
	c3TitleRGB = getTitle();
	run("RGB Color");
	Ext.Manager3D_MonoSelect();
	Ext.Manager3D_DeselectAll();
	for(i=0;i<nb;i++) {
	    Ext.Manager3D_Measure3D(i,"Vol",V);
	    mfr=(V-min)/(max-min);
	    Ext.Manager3D_Select(i);
    	Ext.Manager3D_FillStack(255*mfr, 255, 255*mfr);
	}
	for(i=0;i<nb;i++) {
		selectWindow(c3Title);
		Ext.Manager3D_Quantif3D(i, "Mean", meanM);
	    if (meanM < minMacroValueRed) {
	    	selectWindow(c3TitleRGB);
	    	Ext.Manager3D_Select(i);
    		Ext.Manager3D_FillStack(255, 0, 0);
	    }
	}
	selectWindow(c3TitleRGB);
	run("Set... ", "zoom=150"); 
	setLocation(100, 220);
	getDimensions(w, h, channels, slices, frames);
	m_slice = slices/2;
	for (i=1; i<=slices; i++) {
		Stack.setPosition(1, i, 1);
		wait(120);
	}
	for (i=slices; i>=1; i--) {
		Stack.setPosition(1, i, 1);
		wait(120);
	}	
	Stack.setPosition(1, m_slice, 1);

	open(inputDir + cNr + "-" + iNr + ".tif");
	title4D = getTitle(); 
	selectWindow(title4D);	
	Overlay.show;
	run("Set... ", "zoom=150"); 
	setLocation(900, 220);
	getDimensions(w, h, channels, slices, frames);
	m_slice = slices/2;
	Stack.setPosition(1, m_slice, 1);
	run("Green"); 
	run("Enhance Contrast", "saturated=0.5");
	Stack.setPosition(2, m_slice, 1);
	run("Red"); 
	run("Enhance Contrast", "saturated=0.5");
	Stack.setPosition(3, m_slice, 1);
	run("Blue"); 
	run("Enhance Contrast", "saturated=0.5");
	run("Make Composite");
	for (i=1; i<=slices; i++) {
		Stack.setPosition(1, i, 1);
		wait(80);
	}
	for (i=slices; i>=1; i--) {
		Stack.setPosition(1, i, 1);
		wait(80);
	}	
	Stack.setPosition(1, m_slice, 1);
	selectWindow(c3TitleRGB);

	if (nb > 5) {
		selectWindow(c3TitleRGB);
		for (i=1; i<=slices; i++) {
			Stack.setPosition(1, i, 1);
			wait(120);
		}
		for (i=slices; i>=1; i--) {
			Stack.setPosition(1, i, 1);
			wait(120);
		}	
		Stack.setPosition(1, m_slice, 1);
	}
	if (repeatcomposite) {
		selectWindow(title4D);
		for (i=1; i<=slices; i++) {
			Stack.setPosition(1, i, 1);
			wait(120);
		}
		for (i=slices; i>=1; i--) {
			Stack.setPosition(1, i, 1);
			wait(120);
		}	
		Stack.setPosition(1, m_slice, 1);
		selectWindow(c3TitleRGB);
	}
	
	run("Brightness/Contrast...");

	call("ij.gui.WaitForUserDialog.setNextLocation",x,y);
	waitForUser("Deselect unwanted objects.");
	
	//Format:	cNr	iNr	label	x	y	z	pSize	pixels	maxDiameter	minDiameter	roundness	lysosome	macrophage
	resultF = outputDir + cNr + "-" + iNr + "-" + "result.txt";
	f = File.open(resultF);
	Ext.Manager3D_Count(nb_obj);
	for (object=0; object < nb_obj; object++) {
		selectWindow(c2Title);
		//The object number and the type of measure ("IntDen", "Mean", "Min", "Max", "Sigma") 
		Ext.Manager3D_Quantif3D(object, "Mean", mean);;
		Ext.Manager3D_Centroid3D(object, cx, cy, cz);
		//The object number and the type of measure ("Vol", "Surf", “NbVox”, "Comp", "Feret", "Elon1", "Elon2", "DCMin", "DCMax", "DCMean", "DCSD")	
		Ext.Manager3D_Measure3D(object,"Vol",volume);
		Ext.Manager3D_Measure3D(object,"NbVox",nbVox);
		Ext.Manager3D_Measure3D(object,"Elon1",elon1);
		Ext.Manager3D_Measure3D(object,"Elon2",elon2);
		selectWindow(c3Title);
		//The object number and the type of measure ("IntDen", "Mean", "Min", "Max", "Sigma") 
		Ext.Manager3D_Quantif3D(object, "Mean", meanM);
		selectWindow(titleSeg);
		Ext.Manager3D_Quantif3D(object, "Mean", meanSeg);
		print(f, cNr + "4\t" + iNr + "\t" + meanSeg + "\t" + (cx*width) + "\t" + (cy*height) + "\t" + (cz*depth) + "\t" + volume + "\t" + nbVox + "\t" + elon1 + "\t" + elon2 + "\t" + 1 + "\t" + d2s(mean,3) + "\t" + meanM + "\n");
	}
	objectsFile = outputDir + cNr + "-" + iNr + "-objects.zip";
	Ext.Manager3D_Save(objectsFile);
	File.close(f);
	Ext.Manager3D_Close();
	
//	run("Sync Windows");
	
	while (nImages>0) { 
		selectImage(nImages); 
		close(); 
	}

	filenumber++;
}


function drawEllipseL(x, y, color, lineWidth) {
  	setColor(color);
    setLineWidth(lineWidth);
    Overlay.drawEllipse(x, y, 60, 60);
	}

function drawEllipseS(x, y, color, lineWidth) {
  	setColor(color);
    setLineWidth(lineWidth);
    Overlay.drawEllipse(x, y, 40, 40);
	}

function drawEllipseM(x, y, color, lineWidth) {
  	setColor(color);
    setLineWidth(lineWidth);
    Overlay.drawEllipse(x, y, 50, 50);
	}
