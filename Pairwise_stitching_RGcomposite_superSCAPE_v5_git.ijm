/* Stitching pipeline for roving scans
 * 
 * optimized for mediSCAPE data
 *  --- Hillmanlab x 2021 ---
 * For continuous stitching of roving SCAPE scans when canvasing a large region of tissue. 
 * This requires:
 * 	- .tiff stacks for each color and timepoint
 * optional: flatfield correction matrices for the data (here "gradC#.tiff" - needs to be opened upon running script)
 * 
 */

setBatchMode(true); // if true, stitching works in the background of ImageJ

// Set up file handling (path and filename) --- NEED to adapt!
fpath = "E:/SCAPE_data/tiffs/";	// exemplatory Windows path to folder of .tiff stacks to stitch
fname = "mouth_lip_run1_t";		// exemplatory file stem (used for all .tiffs to stich in this folder & running over time t)

// Set up stitching paramenter --- adapt to dataset characteristics (fast rove -> small nvol; large FOV -> bigger nvol)
start = 1; 		// starting stitching at volume t = start
endt = 360; 	// end stitching at volume t = endt
nvol = 1;  		// how many volumes to skip in between stitched volumes. If the R value is not above Rmin, consecutive ("attempt") volumes will be tried until R>Rmin
Rmin = 0.90; 	// minimum correlation coefficient (goodness of stitch) for a successful stitch
Rfloor = 0.70; 	// if Rmin was not reached, consecutive vollumes will be attempted and their highest R-value will be stitched if R>Rfloor
attempts = 5; 	// number of attempts until either one with R>Rmin is found or the best if R>Rfloor


// Initialize
num = nvol+start;	// num is the timepoint t of the next stitching candidate .tiff-stack
numold = start;		// set the base volume to start stitching to
count = 0;			// counts stitching attempts based on current intermediate result
count2 = 0;			// counts console outputs (used for clearing them)
lasttif = fpath + "Fused_" + numold + ".tif"; // if stitching should be continued from a previous stitching run (this is optional and would be dialed in later)


// Initialize positions for volumes to initially stitch at
x = 0;
y = 0;
z = 0;
x1 = x;
y1 = y;
x_prev = 0;
y_prev = 0;
z_prev = 0;
dx = 0;		// delta for initializing registration shifts by previous run
dy = 0;
dz = 0;
er1 = 0; 	// permitted decrease of Rmin after unsuccessful stitching trys

// miniSCAPE specifics: crop the Y and X for stitching (small and round FOV of miniSCAPE)
//maxw = 450;
//maxh = 110;//90;


// --- load initial volumes and merge them to RG volumes ---
// 		filenaming convention "filestem"_t#_C#.tiff for timepoinnt t and color-channel C
// 		NEED TO ADAPT: naming convention and crop parameter
// load channel 1 (red)
open(fpath+fname+IJ.pad(start, 7)+"_C1.tiff");											// load .tiff stack from file location
getDimensions(w0, h0, channels, slices, frames);										// get the dimensions (for cropping boundaries)
run("Slice Remover", "first="+(slices-4)+" last="+slices+" increment=1"); wait(100);	// remove last 4 slices
run("Slice Remover", "first=1 last=50 increment=1");									// remove first 50 slices
rename("AR"); wait(100);																// rename the Fiji window (important for window handling)
//imageCalculator("Multiply stack", "AR","gradC1.tiff"); wait(100);						// optional flatfield corrention using gradient image "gradC1.tiff"

// load channel 2 (green)
open(fpath+fname+IJ.pad(start, 7)+"_C2.tiff");
run("Slice Remover", "first="+(slices-4)+" last="+slices+" increment=1"); wait(100);	// remove last 4 slices
run("Slice Remover", "first=1 last=50 increment=1");									// remove first 50 slices
rename("AG"); wait(100);
//imageCalculator("Multiply stack", "AG","gradC2.tiff"); wait(100);

// merge both channels (Composite 1) as the volume to stitch to
run("Merge Channels...", "c1=AR c2=AG create"); wait(100);
//makeRectangle(round((w0-maxw)/2), 10, maxw, maxh);
//run("Crop"); wait(100);
run("Subtract...", "value=105 stack"); // background/offset removal
makeRectangle(0, 0, w0, h0); 
rename("Composite1"); wait(100);
getDimensions(w, h, channels1, slices1, frames1);

//makeRectangle(round((w-maxw)/2), round((h-maxh)/2), maxw, maxh); // define ROI for stitching


// --- enter into for-loop for pairwise stitching ---
while (num<(endt+1)) {

	// reset R values for this stitching run
	R = 0;				 // current stitching attempt
	R_bestAttempt = 0;	 // in case further volumes need to be tried for stitching save R and num of the best attempt
	num_bestAttempt = 0; // see above
	
	// --- load next stack ( dual color ) to try stitching (alternatively, load from "lasttif")
	print("  ----- Trying to stitch to the current volume "+numold+" the next volume "+num+"  -----");
	// load next C2 (green) stack
	open(fpath+fname+IJ.pad(num, 7)+"_C2.tiff");
	run("Slice Remover", "first="+(slices-4)+" last="+slices+" increment=1"); wait(100);	// remove last 4 slices
	run("Slice Remover", "first=1 last=50 increment=1");									// remove first 50 slices
	rename("AG"); wait(100);
	//imageCalculator("Multiply stack", "AG","gradC2.tiff"); wait(100);

	// load next C1 (red) stack
	open(fpath+fname+IJ.pad(num, 7)+"_C1.tiff");
	run("Slice Remover", "first="+(slices-4)+" last="+slices+" increment=1"); wait(100);	// remove last 4 slices
	run("Slice Remover", "first=1 last=50 increment=1");									// remove first 50 slices
	rename("AR"); wait(100);
	//imageCalculator("Multiply stack", "AR","gradC1.tiff"); wait(100);

	// merge channels for dual color (Composite2)
	run("Merge Channels...", "c1=AR c2=AG create"); wait(100);
	//makeRectangle(round((w0-maxw)/2), 10, maxw, maxh); 
	//run("Crop"); wait(100);
	run("Subtract...", "value=105 stack");
	makeRectangle(0, 0, w0, h0); 
	rename("Composite2"); wait(100);
	
	//makeRectangle(round((w-maxw)/2), round((h-maxh)/2), maxw, maxh); // define ROI for stitching

	// try the stitching and obtain possible R value (cast from text)
	run("Pairwise stitching", "first_image=Composite1 second_image=Composite2 fusion_method=[Do not fuse images] fused_image=[Composite3] check_peaks=50 ignore compute_overlap x="+xd+" y="+yd+" z="+zd+" registration_channel_image_1=[Average all channels] registration_channel_image_2=[Average all channels]"); wait(500);

	logString = getInfo("log");	// read the R value from console log
	rIndex1 = lastIndexOf(logString, "(R)=");
	R = substring(logString, rIndex1+4, rIndex1+9);


	// stitch if good enough
	if(R > Rmin-er1){

		// save previous shifts
		x_prev = x; y_prev = y; z_prev = z;
		
		// cast stitching parameter (next stack relative xyz-positioning)
  	 	xIndex1 = lastIndexOf(logString, "first): (");
	  	yIndex1 = indexOf(substring(logString, xIndex1+9,xIndex1+40), ", ")+xIndex1+9;
	  	zIndex1 = indexOf(substring(logString,yIndex1+2, yIndex1+40), ", ")+yIndex1+2;
	  	endIndex = indexOf(substring(logString,zIndex1+2, zIndex1+40), ")")+zIndex1+2;

	  	x = parseFloat(substring(logString, xIndex1+9,yIndex1));
	  	y = parseFloat(substring(logString,yIndex1+2, zIndex1));
	  	z = parseFloat(substring(logString,zIndex1+2, endIndex));

		// stitch together (-> Composite3)
		run("Pairwise stitching", "first_image=Composite1 second_image=Composite2 fusion_method=[Linear Blending] fused_image=[Composite3] check_peaks=50 ignore x="+x+" y="+y+" z="+z+" registration_channel_image_1=[Average all channels] registration_channel_image_2=[Average all channels]"); wait(500);

		// save previous stitched volume and use stitching result as new base
		selectWindow("Composite1");
		saveAs("Tiff", fpath+"Fused_"+start+"-"+num); wait(500); close();
		lasttif = fpath+"Fused_"+start+"-"+num+".tif";

		selectWindow("Composite3");
		rename("Composite1");
		getDimensions(wn, hn, channelsn, slicesn, framesn);

		// manage boundary condition for the latest fusing result (wn: width new, hn: height new)
		// and define new ROI
		if ((x1+w)>wn){	x1 = wn-w; };
		if ((y1+h)>hn){	y1 = hn-h; };

		selectWindow("Composite1");
		makeRectangle(round(x1), y1, round(w), h);


	  	// update counts
	  	count = 0;
		numolder = numold;
		numold = num;
		num = num + nvol;
		er1 = 0;
		R_bestAttempt = 0;
		num_bestAttempt = 0;

		// compute the deltas of the previous shifts
		xd = x-x_prev; yd = y-y_prev; zd = 0; //instead of z-z_prev;
  	  	
  	  	// select roi coordinates for stitching next image - if neg set to 0
	  	if (abs(x)==x){	x1 = x;	} else { x1 = 0; };
	    if (abs(y)==y){	y1 = y;	} else { y1 = 0; };
	  	
	  	print("SUCCESS: images fused with R = " +R);
	  	
	}else{ 	// not enough -> don't stitch
      	print("FAIL: image not fused");

		// check if this was so far the best attempt in terms of R-value
		if (R > R_bestAttempt){
			R_bestAttempt = R;
			num_bestAttempt = num;
		};
		print("Best fit so far -> volume " + num_bestAttempt + " with R = " + R_bestAttempt);
		
		num = num + 1; //try next consecutive volume
      	count = count +1;
		
		// if all attempts are used up, choose the best try if that was reasonable good
		if (count > attempts){
			if(Rfloor < R_bestAttempt){
				print("-- back to the best fit so far -> vol " + num_bestAttempt + " and R = " + R_bestAttempt);
				num = num_bestAttempt;	// little workaround to set the best try as the one to stitch
				er1 = Rmin - R_bestAttempt + 0.1;	// set the er1, that subtracts from R per round, artificially high to have the next one pass
			}else{
				// stop here and suggest to continue with next result
				print("-- stopped at volume " + numold);
				num = endt+1;	// break the while loop
			};
		};
	}; 		// end of if R > Rmin clause


	// close the volume that was considered to stitch
	selectWindow("Composite2");	close();

	// count console outputs and clear if necessary
	count2 = count2 + 1;
	if (count2>40){
		call("java.lang.System.gc");
		count2 = 0;
		print("\\Clear");
	};


}; // end of while-loop

