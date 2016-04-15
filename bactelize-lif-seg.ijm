
//E data: e_PC
//E length: Sphere_Regularization
//Initializatioin: LocalMax
//Region Competition Parameters
//no fusion, no fission, handles yes
//Lambda E length: 0.0400
//Thete E merge: 0.000
//Max Iterations: 300
//Oscillation threshold (convergence): 0.0200
//Curvature based gradient flow option: 8
//Local Max Initialization:
//Radius: 10
//Sigma: 3
//Tolerance: 0.4
//Region Tol: 20

setBatchMode(true);
x=100; y=130;
call("ij.gui.WaitForUserDialog.setNextLocation",x,y);
waitForUser("Run Region competition before this Macro to set all parameters correct");


//  inputDir = getDirectory("Choose the Directory of the input files. ");
//inputDir = "/media/mmeuli/WD-HD-ext4/20160408_BCG_Pasteur-Aeras_in_THP-1/data-deconvolved/in/test/"
inputDir = "/media/mmeuli/WD-HD-ext4/20160408_BCG_Pasteur-Aeras_in_THP-1/data-deconvolved/in/smeg/"
//  inputDir = getDirectory("Choose the Directory for the output files. ");
//outputDir = "/media/mmeuli/WD-HD-ext4/20160408_BCG_Pasteur-Aeras_in_THP-1/data-deconvolved/out/test/"
outputDir = "/media/mmeuli/WD-HD-ext4/20160408_BCG_Pasteur-Aeras_in_THP-1/data-deconvolved/out/"
list = getFileList(inputDir);
count = 0;
for (i=0; i<list.length; i++) {
	if (endsWith(list[i], "-bacteria.tif")) {
		count++;
	}
}

print("Number of -bacteria.tif files = " + count);
bacFiles = newArray(count);
bacFilesCount = 0;
for (i=0; i<list.length; i++) {
	if (endsWith(list[i], "-bacteria.tif")) {
		bacFiles[bacFilesCount] = inputDir + list[i];
		print((bacFilesCount) + ": " + list[i]);
		bacFilesCount++;
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

while (filenumber < bacFilesCount) {	
	filename = bacFiles[filenumber];
	fileT = File.getName(filename); 
	columns=split(fileT,"-"); 
	cNr=parseInt(columns[0]); 
	iNr = parseInt(columns[1]); 
	print("Opening file: " + filename);
	open(filename);
	title = getTitle(); 
	selectWindow(title);
	run("Region Competition", "e_data=e_PC e_length=Sphere_Regularization lambda=0.0400 theta=0.0000 max_iterations=300 oscillation=0.0200 initialization=LocalMax inputimage=[title] labelimage=[]");
  	titleSeg = getTitle(); 
	selectWindow(titleSeg);
	resultS = outputDir + cNr + "-" + iNr + "-" + "labeled.tif";
	saveAs("Tiff", resultS);
	close();
	selectWindow(title);
	close();
	filenumber++;
}





