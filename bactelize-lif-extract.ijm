setBatchMode(true);
path = File.openDialog("Select a File");
x=100; y=130;
Dialog.create("Prefix");
Dialog.addNumber("Enter prefix (e.g. for wt: 1 :", 1);
Dialog.setLocation(x,y);
Dialog.show();
cNr = Dialog.getNumber();

run("Bio-Formats Macro Extensions");
Ext.setId(path);
Ext.getCurrentFile(file);
Ext.getSeriesCount(seriesCount);
print(seriesCount);


for (s=1; s<=seriesCount; s++) {
	run("Bio-Formats Importer", "open=&path autoscale color_mode=Default view=Hyperstack stack_order=XYCZT series_"+s);  
	inputDir = "/media/mmeuli/WD-HD-ext4/20150325_BCG_Pasteur-Aeras_zmp1_ko_in_RAW/2016014/in/";
	save4D_path = inputDir + cNr + "-" + s + ".tif";   
	saveAs("tiff", save4D_path);
	title = getTitle();
	run("Split Channels");
	c1Title = "C1-" + title; 
	c2Title = "C2-" + title; 
	c3Title = "C3-" + title; 
	selectWindow(c1Title);
	resultB = inputDir + cNr + "-" + s + "-" + "bacteria.tif";
	saveAs("Tiff", resultB);
	selectWindow(c2Title);
	resultL = inputDir + cNr + "-" + s + "-" + "lysosomes.tif";
	saveAs("Tiff", resultL);
	selectWindow(c3Title);
	resultM = inputDir + cNr + "-" + s + "-" + "macrophage.tif";
	saveAs("Tiff", resultM);;
}



//run("Make Composite");
//run("Brightness/Contrast...");


//out_path = getDirectory("image") + getTitle() + ".tif";