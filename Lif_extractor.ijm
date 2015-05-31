
macro "LIF Extractor" {


print("\\Clear");
//print("DIR_PATH :"+DIR_PATH);






// Creation of the dialog box
Dialog.create("LIF Ectractor");
Dialog.addMessage("Select folder containing lif files.\n");
Dialog.addMessage("Macro extracts lif files as tiff stacks.\n");
 

// Feeding variables from dialog choices

SPLIT=0;


// Get the folder name 
DIR_PATH=getDirectory("Select a directory");


// Get all file names
// or change macro for single file using
//File.openDialog(path)

ALL_NAMES=getFileList(DIR_PATH);
ALL_EXT=newArray(ALL_NAMES.length);
// Create extensions array
for (i=0; i<ALL_NAMES.length; i++) {
	LENGTH=lengthOf(ALL_NAMES[i]);
	ALL_EXT[i]=substring(ALL_NAMES[i],LENGTH-4,LENGTH);
}


setBatchMode(true);




// Loop on all .lei and .lif extensions
for (n=0; n<ALL_EXT.length; n++) {
if (ALL_EXT[n]==".lei" || ALL_EXT[n]==".lif") {

	
	// Get the file path
	FILE_PATH=DIR_PATH+ALL_NAMES[n];
	
	// Store components of the file name
	FILE_NAME=File.getName(FILE_PATH);
	FILE_PATH_LENGTH=lengthOf(FILE_PATH);
	FILE_NAME_LENGTH=lengthOf(FILE_NAME);
	FILE_DIR=substring(FILE_PATH,0,FILE_PATH_LENGTH-FILE_NAME_LENGTH);
	FILE_EXT=substring(FILE_NAME,FILE_NAME_LENGTH-4,FILE_NAME_LENGTH);
	FILE_SHORTNAME=substring(FILE_NAME,0,FILE_NAME_LENGTH-4);

print("Extracting file :", FILE_NAME);	

//print("");	
//print("FILE_PATH:", FILE_PATH);
//print("FILE_NAME:", FILE_NAME);	
//print("FILE_DIR:", FILE_DIR);
//print("FILE_EXT:", FILE_EXT);
//print("FILE_SHORTNAME:", FILE_SHORTNAME);

	
	// Localize or create the output folder
	OUTPUT_DIR="Void";
	OUTPUT_DIR=FILE_DIR+FILE_SHORTNAME+"_ZStacks"+File.separator;
	File.makeDirectory(OUTPUT_DIR);
	


	
//print("OUTPUT_DIR: "+OUTPUT_DIR);
//print("");
	
	
	// Start BioFormats and get series number in file.
	run("Bio-Formats Macro Extensions");
	Ext.setId(FILE_PATH);
	Ext.getSeriesCount(SERIES_COUNT);
	SERIES_NAMES=newArray(SERIES_COUNT);

print(SERIES_COUNT + " images found");	
//print(SERIES_COUNT + " images found in "+ ALL_EXT[n] +" file");
print("");	
	// Loop on all series in the file
	for (i=0; i<SERIES_COUNT; i++) {
		
		// Get serie name and channels count
		Ext.setSeries(i);
		Ext.getEffectiveSizeC(CHANNEL_COUNT);
		SERIES_NAMES[i]="";
		Ext.getSeriesName(SERIES_NAMES[i]);
		TEMP_NAME=toLowerCase(SERIES_NAMES[i]);
		
//print("SERIES_NAMES["+i+"]: "+ SERIES_NAMES[i] + " (TEMP_NAME: " + TEMP_NAME +")");

		// Import the serie (split channels)
//		run("Bio-Formats Importer", "open=["+ FILE_PATH + "] " + "view=[Standard ImageJ]" + " stack_order=Default split_channels " + TEMP_NAME);
//		print("Bio-Formats Importer", "open=["+ FILE_PATH + "] " + "split_channels view=[Standard ImageJ]" + " stack_order=Default " + "series_"+d2s(i+1,0));
		
		//run("Bio-Formats Importer", "open=["+ FILE_PATH + "] " + "split_channels view=[Standard ImageJ]" + " stack_order=Default " + "series_"+d2s(i+1,0));
	if(!SPLIT){
		run("Bio-Formats Importer", "open=["+ FILE_PATH + "] " + "autoscale view=[Hyperstack]" + " stack_order=Default " + "series_"+d2s(i+1,0));
	}
	else{
	run("Bio-Formats Importer", "open=["+ FILE_PATH + "] " + "split_channels view=[Standard ImageJ]" + " stack_order=Default " + "series_"+d2s(i+1,0));
	}
			
			// Construct window name
			//TEMP_CHANNEL=d2s(j,0);
			// Windows has Series Name in title only if more than one Serie
			if(!SPLIT){
			if(SERIES_COUNT==1) {
				SOURCE_WINDOW_NAME=FILE_NAME;//+ " - C="+TEMP_CHANNEL;
			}
			else {
			SOURCE_WINDOW_NAME=FILE_NAME+" - "+SERIES_NAMES[i];//+" - C="+TEMP_CHANNEL;
			}
			TYPE="";
			
			//Select source image and filter if asked
			selectWindow(SOURCE_WINDOW_NAME);
			
			
			// rename image according to processing
			NEW_WINDOW_NAME=FILE_NAME+" - "+SERIES_NAMES[i];
			rename(NEW_WINDOW_NAME);

//print("NEW_WINDOW_NAME: "+NEW_WINDOW_NAME);

		
			// Create output file path and save the output image
			OUTPUT_PATH="Void";
				NEW_WINDOW_NAME=replace(NEW_WINDOW_NAME, '/', '_');
print("Image " + i+1 + ":" + NEW_WINDOW_NAME);
				OUTPUT_PATH=OUTPUT_DIR+NEW_WINDOW_NAME+".tif";
				save(OUTPUT_PATH);
				close();
			}
			else{
				// Loop on each channel (each opened window)
				for(j=0; j<CHANNEL_COUNT; j++) {
					TEMP_CHANNEL=d2s(j,0);
					if(SERIES_COUNT==1) {
						SOURCE_WINDOW_NAME=FILE_NAME+ " - C="+TEMP_CHANNEL;
					}
					else {
					SOURCE_WINDOW_NAME=FILE_NAME+" - "+SERIES_NAMES[i]+" - C="+TEMP_CHANNEL;
					}
					TYPE="";
				//Select source image and filter if asked
				selectWindow(SOURCE_WINDOW_NAME);
			
			
				// rename image according to processing
				NEW_WINDOW_NAME=FILE_NAME+" - "+SERIES_NAMES[i]+" - C="+TEMP_CHANNEL;
				rename(NEW_WINDOW_NAME);

//print("NEW_WINDOW_NAME: "+NEW_WINDOW_NAME);

		
				// Create output file path and save the output image
				OUTPUT_PATH="Void";
				NEW_WINDOW_NAME=replace(NEW_WINDOW_NAME, '/', '_');
print("Image " + i+1 + ":" + NEW_WINDOW_NAME);				
				OUTPUT_PATH=OUTPUT_DIR+NEW_WINDOW_NAME+".tif";
				save(OUTPUT_PATH);
				close();
				
				}//end channel loop
				
			}//end else statement : channel splitting

			
// print("OUTPUT_PATH :"+OUTPUT_PATH);
			


			
				
	}
print("");
}	// end of IF test on lei and lif extensions
}	// end of FOR loop on n extensions
print("Done");
setBatchMode("exit and display");
showStatus("finished");
}	// end of macro
	