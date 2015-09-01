macro "Open RC-bactelize-Out-Images" {

pixelW=0.0616
pixelH=0.0616
pixelD=0.1998
input = getDirectory("Directory containing images to open");

Dialog.create("File prefix");
Dialog.addString("File prefix: ", "1-1", 5);
Dialog.show();
prefix = Dialog.getString();

pathfileSub2="";
pathfileSub2=File.openDialog("Choose the dataset for selection of bacteria:"); 
//pathfileSub2="/home/mmeuli/batch/out/1-1-A_result-subset2.txt";
//pathfileSub2="/home/mmeuli/Bioimage/Colocalization-Experiments/20150611 BCG Pasteur-Aeras zmp1 ko in RAW/data-results-deconvoluted/20150818-results/1-1-A_result-subset2.txt";

cNr=0;
iNr=0;
columns=split(prefix,"-"); 
cNr=parseInt(columns[0]); 
iNr=parseInt(columns[1]); 

while (nImages>0) { 
	selectImage(nImages); 
	close(); 
	} 

open(input + prefix + "-bacteria.tif");
run("Properties...", "pixel_width=pixelW pixel_height=pixelH voxel_depth=pixelD");
run("Clear Results");
run("Statistics"); 
min=getResult("Min", 0);
max=getResult("Max", 0);
setMinAndMax(min, max);
run("Set... ", "zoom=60 x=416 y=416");
setLocation(1300, 25);

open(input + prefix + "-lysosomes.tif");
run("Properties...", "pixel_width=pixelW pixel_height=pixelH voxel_depth=pixelD");
run("Clear Results");
run("Statistics"); 
min=getResult("Min", 0);
max=getResult("Max", 0);
setMinAndMax(min, max);
run("Set... ", "zoom=60 x=416 y=416");
setLocation(740, 620);

open(input + prefix + "-macrophage.tif");
run("Properties...", "pixel_width=pixelW pixel_height=pixelH voxel_depth=pixelD");
run("Clear Results");
run("Statistics"); 
min=getResult("Min", 0);
max=getResult("Max", 0);
setMinAndMax(min, max);
run("Set... ", "zoom=60 x=416 y=416");
setLocation(1300, 620);

open(input + prefix + "-labeled.tif");
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

  
function drawEllipse(x, y, color, lineWidth) {
  	setColor(color);
    setLineWidth(lineWidth);
    Overlay.drawEllipse(x, y, 100, 100);
	}


}	// end of macro