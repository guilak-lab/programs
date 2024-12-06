//select folder where .czi files are and folder where .csv result files will be saved
input = getDirectory("Choose input directory"); 
output = getDirectory("Choose output directory"); 

setBatchMode(true);
//name a variable to hold the list of the input files
list = getFileList(input);

//for-loop for iterating through the input list of files and running the analysis
for(i=0;i < list.length; i++){
	if (endsWith(list[i], ".czi")) {
		open(input + list[i]);
		MeIF (input, output, list[i]);
	}
}
setBatchMode(false);

//open file as composite and do not split channels		
function MeIF (input, output, filename) {

//split the channels for file 
	run("Split Channels");
	file_title = File.getNameWithoutExtension(list[i]);

//new channel file for 488 (H3K9Me3 IF), generate MIP 2D image and rename it as the name of the channel
	selectWindow("C1-" + file_title + ".czi"); 
	run("Z Project...", "projection=[Max Intensity]"); 
	run("Median...", "radius=10");
	rename("channel-488");
		
	selectWindow("C2-" + file_title + ".czi"); 
	run("Z Project..." , " projection = [Max Intensity]");
	run("Median...", "radius=10");
	rename("channel-405");

//create the duplicate of the 405 MIP to use as mask for ROI

	run("Duplicate...", "title=DAPI_ROI");

//auto thresholds the DAPI_RIO image; manually you select Image>Adjust>AutoThreshold> Triangle Method
	selectWindow("DAPI_ROI"); 
	setAutoThreshold("IsoData dark no-reset");
//setThreshold(54, 255);
	setOption("BlackBackground", true);
//converts it to binary black0-white255(BW) mask; Process>Binary>Convert to mask
	run("Convert to Mask");
	
//removes background pixels, splits cells, fills holes within cells
	run("Remove Outliers...", "radius=2 threshold=50 which=Bright");
	run("Watershed");
		
//select the measurements you want to measure (mean = Mean Gray Value); note "display" will display the name of channel
	run("Set Measurements...", "mean display decimal=3");

//analyze the ROI image
	selectWindow("DAPI_ROI");
	run("Analyze Particles...", "size=500-Infinity pixel circularity=0.75-1.00 display exclude clear add in_situ");
	close("Results");

//analyze 488 channel
	selectWindow("channel-488");
	roiManager("Show None");
	roiManager("Show All");
	roiManager("Measure");

//save results as CSV file
	saveAs("Results", output + file_title + "_results.csv");
	run("Close All");
		
	}