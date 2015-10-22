macro "Open RC-bactelize-Out-Images" {

pixelW=0.0499610
pixelH=0.0499610
pixelD=0.1502000

path-to-images = getDirectory("Directory containing images to open");

Dialog.create("File prefix");
Dialog.addString("File prefix: ", "1-1", 5);
Dialog.show();
prefix = Dialog.getString();

path-to-dataset="";
path-to-dataset=File.openDialog("Choose the dataset for selection of bacteria:"); 
//path-to-dataset="/home/mmeuli/batch/out/1-1-A_result-subset2.txt";
//path-to-dataset="/home/mmeuli/Bioimage/Colocalization-Experiments/20150611 BCG Pasteur-Aeras zmp1 ko in RAW/data-results-deconvoluted/20150818-results/1-1-A_result-subset2.txt";

cNr=0;
iNr=0;
columns=split(prefix,"-"); 
cNr=parseInt(columns[0]); 
iNr=parseInt(columns[1]); 

while (nImages>0) { 
	selectImage(nImages); 
	close(); 
	} 

open(path-to-images + prefix + "-bacteria.tif");
run("Properties...", "pixel_width=pixelW pixel_height=pixelH voxel_depth=pixelD");
run("Clear Results");
run("Statistics"); 
min=getResult("Min", 0);
max=getResult("Max", 0);
setMinAndMax(min, max);
run("Set... ", "zoom=60 x=416 y=416");
setLocation(1300, 25);

open(path-to-images + prefix + "-lysosomes.tif");
run("Properties...", "pixel_width=pixelW pixel_height=pixelH voxel_depth=pixelD");
run("Clear Results");
run("Statistics"); 
min=getResult("Min", 0);
max=getResult("Max", 0);
setMinAndMax(min, max);
run("Set... ", "zoom=60 x=416 y=416");
setLocation(740, 620);

open(path-to-images + prefix + "-macrophage.tif");
run("Properties...", "pixel_width=pixelW pixel_height=pixelH voxel_depth=pixelD");
run("Clear Results");
run("Statistics"); 
min=getResult("Min", 0);
max=getResult("Max", 0);
setMinAndMax(min, max);
run("Set... ", "zoom=60 x=416 y=416");
setLocation(1300, 620);

open(path-to-images + prefix + "-labeled.tif");
run("Properties...", "pixel_width=pixelW pixel_height=pixelH voxel_depth=pixelD");
run("Clear Results");
run("Statistics"); 
min=getResult("Min", 0);
max=getResult("Max", 0);
setMinAndMax(min, max);
run("3-3-2 RGB");
run("Set... ", "zoom=60 x=416 y=416");
setLocation(740, 25);

run("Synchronize Windows");


if (path-to-dataset != ""){
	filestring=File.openAsString(path-to-dataset); 
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

  
function drawEllipse(x, y, color, lineWidth) {
  	setColor(color);
    setLineWidth(lineWidth);
    Overlay.drawEllipse(x, y, 100, 100);
	}


}	// end of macro
