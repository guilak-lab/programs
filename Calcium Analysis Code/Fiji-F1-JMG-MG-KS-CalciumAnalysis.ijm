//Code for Fluo4/FuraRed ratiometric analysis; generates spreadsheet where #rows=#frames and #columns=#cells detected by ROI based on FuraRed channel MIP 
//select folder where .czi files are and folder where .csv result files will be saved
input = getDirectory("Choose input directory"); 
output = getDirectory("Choose output directory"); 

setBatchMode(true); //change to "false" to see windows as macro is running
//name a variable to hold the list of the input files
list = getFileList(input);

//for loop for iterating through the input list of files and running the analysis, function called CaF
for(i=0;i < list.length; i++){
	if (endsWith(list[i], ".czi")) {
		open(input + list[i]);
		CaF (input, output, list[i]);
	}
}
setBatchMode(false);

//on dialog window, open file as "Composite" and do not split channels		
function CaF (input, output, filename) {

//split the channels for file 
	run("Split Channels");
	file_title = File.getNameWithoutExtension(list[i]);

//duplicate, apply median filter and rename channels
	selectWindow("C1-" + file_title +".czi");
	run("Duplicate...", "duplicate");
	run("Median...", "radius=2 stack");
	FuraRed = getTitle();

	selectWindow("C2-" + file_title +".czi");
	run("Duplicate...", "duplicate");
	run("Median...", "radius=2 stack");
	Fluo4 = getTitle();
	
//create another duplicate with median filter of the FuraRed channel for ROI mask	
	selectWindow("C1-" + file_title +".czi");
	run("Duplicate...", "duplicate");
	run("Median...", "radius=2 stack");

//create a maximum intensity projection across stack (i.e. all frames) for newly duped FuraRed channel and rename
	run("Z Project...", "projection=[Max Intensity]");
	FuraRedROI = getTitle();

//Make ROI mask based on FuraRedROI and measure it to add it to the ROI Manager
//auto thresholds the FuraRedROI image; manually you select Image>Adjust>AutoThreshold> Triangle Method
	run("Auto Threshold", "method=IsoData white");
//Watershed binary to separate cells close together
	run("Watershed");
//measure the ROI image to generate mask and add it to the ROI manager,
	run("Analyze Particles...", "size=80-Infinity pixel circularity=0.80-1.00 show=Outlines display exclude clear add in_situ slice");
	close("Results");

//calculate the ratio of Fluo4:FuraRed by dividing the Fluo4/FuraRed
	imageCalculator("Divide create 32-bit stack", Fluo4, FuraRed);
//select the measurements you want to measure and display in results window (i.e. Mean Gray Value); note "display" will display the name of channel
	run("Set Measurements...", "mean decimal=3");
//measure particles in result window of Fluo4/FuraRed ("Result of Fluo4") based on FuraRedROI
//project ROI mask on to result window of Fluo4/FuraRed	
	roiManager("Show None");
	roiManager("Show All");

//measure Mean Gray Value within ROIs across all frames (i.e. through the stack); generate .csv file with columns corresponding to the number of cells detected by ROI and row corresponding to the number of frames (e.g. 300 cells = 300 columns, 140 framed = 140 rows)
	roiManager("Multi Measure measure all");

//save results as a CSV file 
	saveAs("Results", output + file_title + "_results.csv");
	run("Close All");
}
