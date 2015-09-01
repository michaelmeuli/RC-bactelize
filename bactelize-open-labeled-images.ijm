macro "Open Labeled-Images" {

pixelW=0.0616
pixelH=0.0616
pixelD=0.1998
input = getDirectory("Directory containing images to open");


Dialog.create("File type");
Dialog.addString("File suffix: ", "labeled.tif", 20);
Dialog.show();
suffix = Dialog.getString();

pathfileSub2="";
pathfileSub2=File.openDialog("Choose the dataset for selection of bacteria:"); 

cNr=0;
iNr=0;

processFolder(input);

function processFolder(input) {
	list = getFileList(input);
	for (i = 0; i < list.length; i++) {
		if(File.isDirectory(list[i]))
			processFolder("" + input + list[i]);
		if(endsWith(list[i], suffix))
			processFile(input, list[i]);
	}
}

function processFile(input, file) {
	columns=split(file,"-"); 
	cNr=parseInt(columns[0]); 
	iNr=parseInt(columns[1]); 
	open(input + file);
	run("Properties...", "pixel_width=pixelW pixel_height=pixelH voxel_depth=pixelD");
	run("Clear Results");
	run("Statistics"); 
	min=getResult("Min", 0);
	max=getResult("Max", 0);
	setMinAndMax(min, max);
	run("3-3-2 RGB");
	run("Set... ", "zoom=100 x=416 y=416");  //TODO: getImageSize? for x and y
	setLocation(740, 25);

	if (pathfileSub2 != ""){
		filestring=File.openAsString(pathfileSub2); 
		rows=split(filestring, "\n");
		count=0;
		for(i=1; i<rows.length; i++){ 
			columns=split(rows[i],"\t"); 
			if ((parseInt(columns[0])==cNr) && (parseInt(columns[1])==iNr)) {
				count++;
				}
		}
	
		x=newArray(count); 
		y=newArray(count); 
		index=0;
		for(i=1; i<rows.length; i++){ 
			columns=split(rows[i],"\t"); 
			if ((parseInt(columns[0])==cNr) && (parseInt(columns[1])==iNr)) {
				x[index]=parseFloat(columns[3]); 
				y[index]=parseFloat(columns[4]); 
				index++;
				}
			}
		
		Overlay.remove;
		for(i=0; i<count; i++){ 
		 	drawEllipse((x[i]/pixelW-50), (y[i]/pixelH-50), "red", 4);
			Overlay.show;
			}
	}

	
}


function drawEllipse(x, y, color, lineWidth) {
  	setColor(color);
    setLineWidth(lineWidth);
    Overlay.drawEllipse(x, y, 100, 100);
	}


}	// end of macro