input = getDirectory("Input directory");
Dialog.create("File type");
Dialog.addString("File prefix: ", "1-1", 5);
Dialog.show();
prefix = Dialog.getString();
open(input + prefix + "-labeled.tif");

run("3-3-2 RGB");
run("Enhance Contrast...", "saturated=0.1 process_all use");


open(input + prefix + "-bacteria.tif");

run("Enhance Contrast...", "saturated=0.1 process_all use");
