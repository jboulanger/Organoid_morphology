#@File(label="Input file") filename
#@Boolean(label="manual segmentation", value=false) manualsegmentation
#@String(label="Overlay", choices={"Curvature","DNE"},style="listBox") overlaytype
#@Boolean(label="Display info", value=true) iwantinfo
#@Boolean(label="Add colorbar", value=true) iwantcolorbar
#@Float(label="Colorbar min", value = -0.05) colorbarmin
#@Float(label="Colorbar max", value = 0.05) colorbarmax
#@Boolean(label="Save image as jpg and close", value=true) finalclear
#@Boolean(label="Save ROI", value=true) saveroi
#@File(label="Output folder",style="directory") output_folder
#@Boolean(label="Used saved ROI", value=false) useSavedROI
#@Float(label="Marker scale", value=4) markerscale

/*
 * Quantify organoids morphology
 * 
 * usage: press run or batch
 * 
 * Set the parameters
 *  - set the filename
 *  - tick manual segmentation if you want to select manually the outline of the orgnanoid
 *  - set a threshold for segmentation (higher stricter)
 *  - set the minimum area to filter out smaller objects
 *  - tick the last option to save a jpg and close the image
 * 
 * Note: 
 * - If a roi file is present the roi will be use 
 * 
 * Output:
 * - Area/size (need to take into account pixel size)
 * - Shape measure (need to implement other methods)
 *   o fractal dimension
 *   o number of inflection points based on curvature
 *   o Dirichlet Normal Energy (DNE)
 * - Transparency measurement based on variations of intensity
 * 
 * https://onlinelibrary.wiley.com/doi/full/10.1002/ajpa.21489
 * https://arxiv.org/abs/1201.3097
 * 
 * Jerome Boulanger for Ilaria (2020)
 * 
 */

function manualBatchProcess() {
	// process a folder containing tif files
	filelist = getFileList(folder);
	N = filelist.length;
	N = 1;
	for (i = 0; i < N; i++) {
		showProgress(i, N);
	    if (endsWith(filelist[i], ".tif")) { 
			processFile(i, folder, filelist[i], "Organoid Morphology");
	    } 
	}
}

function initTable(tblname) {
	// Initialize a result table if needed and return the next row index.
	if (!isOpen(tblname)) {
		Table.create(tblname);
	}
	selectWindow(tblname);
	return Table.size;
}

function processFile(idx, folder, fname, tblname) {
	// open and process file
	
	setBatchMode("hide");	
	
	roifilename = folder + File.separator + replace(fname, ".tif", ".roi");
	
	open(folder + File.separator + fname);
	id = getImageID();
	
	getPixelSize(unit, pixel_size, pixel_size);
	
	if (useSavedROI && File.exists(roifilename)) {
		roiManager("Open", roifilename);
		roiManager("select", 0);
		fixContourOrientation();
		roiManager("update");
	} else {	
		if (manualsegmentation) {
			if (File.exists(roifilename)) {
				//setBatchMode("exit and display");
				print(roifilename);
				roiManager("Open", roifilename);
				roiManager("select",0);
				run("Interpolate", "interval=50 smooth adjust");	
				run("Fit Spline");
				setTool("freehand");
				setBatchMode("exit and display");
				waitForUser("Please adjust the outline of the organoid and update the roi manager");
				setBatchMode("hide");				
				roiFftSmooth(1);
				run("Fit Spline");
				run("Interpolate", "interval=10 smooth adjust");					
				fixContourOrientation();
				roiManager("add");
				run("32-bit");
			} else {
				setTool("freehand");
				setBatchMode("exit and display");
				waitForUser("Please draw the outline of the organoid");
				setBatchMode("hide");				
				roiFftSmooth(1);
				run("Fit Spline");
				run("Interpolate", "interval=10 smooth adjust");	
				fixContourOrientation();
				roiManager("add");
				run("32-bit");
			}
		} else {						
			segment(id);
		}
	}
	if (roiManager("count")==0) {
		Table.set("File", idx, fname);		
 		print("************"); 		
 		print("** FAILED **");
 		print("************");
		return 0;
	}
	
	// Compute curvature
	roiManager("select", 0);		
	Roi.getCoordinates(x, y);
	run("Select None");
	
	
		
	// compute the length along the contour in micron
	length = getCurvilinearLength(x,y,pixel_size);
	fixNaN(length);
	
	// average radius (distance from contour to barycenter) in micron
	R0 = getAverageRadius(x,y,pixel_size);
	
	// curvature in micron
	curvature = getContourCurvatureGeo(x,y,pixel_size);	
	fixNaN(curvature);
	
	// compute the log of the sum [kappa * length] over the all contour
	E = log(getWeightedEnergy(curvature,length));

	// Detect inflection points such that the curvature is more than 1/(2 R0)
	// or Rc < 2 R0
	inflex = getContourInflectionPoints(curvature, R0/2);
	K = inflex.length;
	
	// Dirichlet normal energy
	dnd = DND(x,y,pixel_size);
	fixNaN(dnd);
	DNE = log(getWeightedEnergy(dnd,length));
	
	// Fractal dimension
	D = roiFractalDimension(x,y,false);
	
	// Transparency
	T = transparency();

	I0 = measureIntensityOutside(id);
		
	//selectWindow(tblname);
	selectImage(id);
	roiManager("select", 0);
	Table.set("File", idx, fname);
 	Table.set("Pixel size ["+unit+"]", idx, pixel_size);
 	
 	Table.set("Area [um^2]", idx, getValue("Area"));
 	Table.set("Perimeter [um]",idx,getValue("Perim."));
 	Table.set("Average radius [um]", idx, R0); 	
 	Table.set("Roundness", idx, getValue("Round"));
 	Table.set("Aspect Ratio", idx, getValue("AR"));
 	Table.set("Feret [um]", idx, getValue("Feret"));
 	Table.set("MinFeret [um]", idx, getValue("MinFeret"));
 	//Table.set("Minor [um]", idx, getValue("Minor"));
 	//Table.set("Major [um]", idx, getValue("Major"));
 	Table.set("Circularity", idx, 100*getValue("Circ."));
 	Table.set("Inflections Points", idx, K);
 	Table.set("Weighted curvature [um^-1]", idx, E);
 	Table.set("DNE", idx, DNE);
 	//Table.set("Fractal Dimension", idx, D);
 	//Table.set("Mean Intensity", idx, getValue("Mean"));
 	//Table.set("StdDev Intensity.", idx, getValue("StdDev"));
 	Table.set("Transparency", idx, T); 	
 	Array.getStatistics(dnd, min, max, mean, stdDev);
 	//Table.set("Min DND", idx, min);
 	//Table.set("Max DND", idx, max);
 	Array.getStatistics(curvature, Cmin, Cmax, Cmean, Cstd);
 	//Table.set("Min curvature [um^-1]", idx, Cmin);
 	//Table.set("Max curvature [um^-1]", idx, Cmax);
 	Table.set("Mean curvature [um^-1]", idx, Cmean);
 	Table.set("Std curvature [um^-1]", idx, Cstd); 	
 	Table.set("R0 x Std curvature", idx,R0*Cstd);
 	//Table.set("Intensity outside", idx, I0);
 	 
 	roiManager("select", 0);
 	//run("Select None");
 	
 	// Overlays
	if (matches(overlaytype,"DNE")) {
		colorContour(x,y,dnd,markerscale);
		if (iwantcolorbar) {
			addColorBar(dnd);
		}
	} else {
		colorContour(x,y,curvature,markerscale);
		if (iwantcolorbar) {
			addColorBar(curvature);
		}
	}
		
	// drawing the inflection points as red circles
	setColor("red");
	for (i = 0; i < K; i++) {
		Overlay.drawEllipse(x[inflex[i]]-markerscale, y[inflex[i]]-markerscale,2*markerscale, 2*markerscale);		
		Overlay.show;
		Overlay.setStrokeWidth(markerscale);
	}

	// Display information as text
	if (iwantinfo) {
		setColor("black");
		setFont("Fixed", 30, "bold");
		Overlay.drawString("TRANSPARENCY: "+d2s(T,3)+"\nR0 * CURVATURE std: "+d2s(R0*Cstd,3)+"\nINFLECTION POINTS: "+K+"\nRADIUS:"+d2s(R0,1)+"um", 0, 30);
		//Overlay.drawString("TRANSPARENCY: "+round(T)+"\nCURVATURE: "+E+"\nDNE: "+DNE+"\nINFLECTION POINTS: "+K+"\nFRACTAL DIM: "+D, 0, 30);
		Overlay.show;
	}
 	 	
 	selectImage(id);
 	run("Select None");
 	if (finalclear) { 
 		ofilename = output_folder + File.separator + replace(fname,".tif",".jpg");
 		print("Saving file \n" + ofilename);
 		saveAs("Jpeg", ofilename);
 		close(); 
 	}
 	
 	if (saveroi) {
 		if (roiManager("count") == 1) {  		
 			ofilename = output_folder + File.separator + replace(fname,".tif",".roi");
	 		print("Saving file \n" + ofilename);
	 		roiManager("select", 0);
	 		roiManager("Save", ofilename);
 		} else {
 			print("There are more than one ROI.. cannot save it as .roi file.");
 		}
 	}
 	 	 	
 	closeRoiManager();
 	run("Select None");
 	setBatchMode("exit and display"); 	
}

function fixContourOrientation() {
	// re-orient the contour so that the normal points outwa	
	Roi.getCoordinates(x, y);
	Array.getStatistics(x, min, max, xc, stdDev);
	Array.getStatistics(y, min, max, yc, stdDev);
	dx = diff(x);
	dy = diff(y);
	dp = 0;
	for (i = 0; i < x.length; i++){
		// dot product of the normal and the radius
		nx = dy[i];
		ny = -dx[i];
		n = sqrt(nx*nx+ny*ny);		
		nx /= n;
		ny /= n;
		rx = x[i]-xc;
		ry = y[i]-yc;
		r = sqrt(rx*rx+ry*ry);
		rx /= r;
		ry /= r;
		if (n > 0 && r > 0) {
			dp += nx*rx+ny*ry;		
		}		
	}	
	dp /= x.length;	
	if (dp > 0) {
		print("Fix contour orientation ("+dp+")");
		Array.reverse(x);
		Array.reverse(y);
		Roi.setPolygonSplineAnchors(x, y);			
	}
}

function fixNaN(array) {
	/* in place replace Not a Number by 0 */
	for (i=0;i<array.length;i++) {		
		if(isNaN(array[i])) {array[i]=0;}
	}
}

function addColorBar(values) {
	// add a colorbar with min/max given by the array values.
	run("RGB Color");
	W = getWidth;
	H = getHeight;
	w = round(0.025*W);
	h = round(0.7*H);
	fs = 20; // font size
	img = getImageID();	
	newImage("Untitled", "8-bit white", w, h, 1);
	cb = getImageID();
	run("Macro...", "code=v=255*(h-y)/h");
	run("Ice");
	run("Select All");
	run("Copy");
	selectImage(img);
	makeRectangle(W-w-3*fs-10, (H-h)/2, w, h);
	run("Paste");
	run("Select All");
	if (colorbarmin == 0 && colorbarmax == 0) {
		Array.getStatistics(values, colorbarmin, colorbarmax, mean, stdDev);
	}
	nticks=5;
	for (i = 0; i < nticks; i++) {
		setFont("Arial", fs);
		setColor("black");
		str = d2s(colorbarmin+i/(nticks-1)*(colorbarmax-colorbarmin), 3);
		Overlay.drawString(str, W-3*fs-5, (H+h)/2-i/(nticks-1)*h);
	}
	Overlay.show();
	selectImage(cb);close();	
}

function segment(id) {
	/* Segment the organoid in the image*/
	removeScaleBar(id);
	preprocessed_id = preprocess2(id);		
	fitContour(preprocessed_id, 1000);	
	selectImage(preprocessed_id); close();	
	//roiManager("select", 0);	
	//roiManager("show none");
	//run("Select None");
}

function removeScaleBar(id) {
	/* Remove the scalebar that was saved onto the image pixels */
	selectImage(id);
	run("Select None");
	a0 = getValue("Area");
	setThreshold(255, 255);	
	run("Create Selection");
	if (getValue("Area")!=a0) {
		if (getValue("Area") != getWidth() * getHeight()) {
			run("Enlarge...", "enlarge=2");
			run("Median...", "radius=10");
			run("Select None");
		}
	}
}

function correctBackground() {
	id0 = getImageID();
	run("Select None");
	run("32-bit");	
	run("Duplicate...", " ");
	id1 = getImageID();
	getDimensions(width, height, channels, slices, frames);
	run("Size...", "width="+width/8+" height="+height/8+" depth=1 constrain average interpolation=Bilinear");
	run("Maximum...", "radius=30");
	run("Gaussian Blur...", "sigma=10");
	run("Size...", "width="+width+" height="+height+" depth=1 constrain average interpolation=Bilinear");
	imageCalculator("subtract", id0, id1);
	selectImage(id1);close();
	selectImage(id0);
}

function preprocess1(id) {	
	/* Provide an initial mask guess */
	print("Preprocessing");
	selectImage(id);
	Overlay.remove();	
	run("Select None");
	getPixelSize(unit, dx, dy);
	run("Duplicate...","title=preprocessed");
	preprocessed_id = getImageID();		
	run("32-bit");
	
	correctBackground();
	for (i = 0; i < 10; i++) {
		run("Minimum...", "radius="+r1);	
		run("Median...",  "radius="+r1);
		run("Maximum...", "radius="+r1);
	}

	selectImage(preprocessed_id);
	run("Duplicate...","title=mask");	
	mask_id = getImageID();
	setAutoThreshold("Li");	
	run("Convert to Mask");		
	r2 = maxOf(1, round(100/dx));	
	run("Minimum...", "radius="+r2);				
	run("Maximum...", "radius="+r2);	
	run("Analyze Particles...", "size=100-Infinity add");
	keepLargestROI();
	roiManager("select",0);
	fixContourOrientation();	
	roiManager("add");
	roiManager("select", 0);
	roiManager("delete");	
	selectImage(mask_id); close();	
	selectImage(preprocessed_id);
	run("Find Edges");
	run("Subtract...", "value=30");
	run("Min...", "value=0");
	return preprocessed_id;
}

function preprocess2(id) {	
	/* Provide an initial mask guess using a contour based approach 
	
	- Add the ROI to the ROI Manager
	- Return the id of the edge image
	*/
	print("Preprocessing");
	selectImage(id);
	closeRoiManager();
	Overlay.remove();	
	run("Select None");
	run("Duplicate...","title=preprocessed");
	preprocessed_id = getImageID();				
	run("32-bit");
	
	getPixelSize(unit, dx, dy);
	r1 = maxOf(1, round(15/dx));		
	for (i = 0; i < 10; i++) {
		run("Maximum...", "radius="+r1);	
		run("Median...", "radius="+r1);
		run("Minimum...", "radius="+r1);
		run("Gaussian Blur...", "sigma=1");	
	}		
	run("Invert");
	idx = getImageID();
	run("Duplicate...","title=edges");
	edges = getImageID();
	run("Find Edges");
	imageCalculator("Multiply", preprocessed_id, edges);
	selectImage(edges); close();
	
	selectImage(preprocessed_id);	
	run("Duplicate...","title=mask");	
	mask_id = getImageID();
	//EMThreshold();
	getStatistics(area, mean, min, max, std, histogram);
	setThreshold(mean+std, max+1);
	//setAutoThreshold("Otsu");
	run("Convert to Mask");		
	
	// close gaps in the contour and fill holes
	getPixelSize(unit, dx, dy);
	r2 = maxOf(1, round(10/dx));
	for (i = 0; i < 10; i++) {
		run("Maximum...", "radius="+r2);	
		run("Median...", "radius="+r2);
		run("Minimum...", "radius="+r2);
	}
	run("Fill Holes");
	
	// cut off small protrusions
	getPixelSize(unit, dx, dy);
	r3 = maxOf(1, round(100/dx));	
	///for (i = 0; i < 5; i++) {
	run("Minimum...", "radius="+r3);
	run("Maximum...", "radius="+r3);	
//	}
	run("Minimum...", "radius="+r3/4);
	run("Analyze Particles...", "size=100-Infinity add");
	keepLargestROI();
	
	roiManager("select",0);
	fixContourOrientation();	
	//run("Interpolate", "interval="+50/dx+"  smooth adjust");
	roiManager("add");
	roiManager("select", 0);
	roiManager("delete");
	
	selectImage(mask_id); close();
	selectImage(preprocessed_id);
	run("Select None");
	getStatistics(area, mean, min, Imax, std, histogram);
	run("Subtract...", "value="+ (Imax / 4));
	run("Min...", "value=0");
	
	return preprocessed_id;
}

function mse(x0,y0,x1,y1) {
	s = 0;
	for (i = 0; i < x0.length; i++) {		
		dmin = 1e6;
		for (j = 0; j < x1.length; j++) {		
			dx = x0[i] - x1[j];
			dy = y0[i] - y1[j];
			d = dx*dx+dy*dy;
			if (d < dmin) {
				dmin = d;
			}
		}
		s += dmin;
	}
	return s / x.length;
}

function fitContour(id, max_iter) {
	/* Fit the contours of the image id, update the ROI manager */		
	print("Fit Contour");	
	
	selectImage(id);
	getPixelSize(unit, dx, dy);	
	
	// Measure the maximum value of the edge image
	run("Select None");
	getStatistics(area, mean, min, Imax, std, histogram);
	
	// Interpolate the contour
	roiManager("select", 0);		
	run("Enlarge...", "enlarge=-2");
	run("Interpolate", "interval="+50/dx+"  smooth adjust");
	
	// Retreive its coordinates and optimize the curve
	Roi.getSplineAnchors(x, y);			
	roiManager("delete");
	
	showStatus("Optimizing contour");
	delta = 1;
	
	//setBatchMode("exit and display");	
	//Roi.setPolygonSplineAnchors(x, y);
	
	x0 = Array.copy(x);
	y0 = Array.copy(y);	
	for (iter = 0; iter < max_iter && delta > 0.1; iter++) {
		showProgress(iter, max_iter);				
		delta = gradTheCurve(x, y, 50/dx, Imax);	
		fftSmooth(x,y,5.);			
		Roi.setPolygonSplineAnchors(x, y);
		run("Interpolate", "interval="+50/dx+" ");
		Roi.getSplineAnchors(x, y);	
		delta = mse(x0,y0,x,y);
		
		x0 = Array.copy(x);
		y0 = Array.copy(y);		
	}
	
	print("Contour fitting finished in "+iter+" iterations.");
	//setBatchMode("exit and display"); exit();
	
	// Update the ROI	
	Roi.setPolygonSplineAnchors(x, y);
	roiManager("add");	
	//exit();
}

function initROI(id) {	
	/* erode the roi and select the biggest resulting region */
	selectImage(id);
	run("Select None");
	run("Duplicate...","title=tmp");	
	run("8-bit");	
	run("Select All");
	setColor(255);
	fill();	
	id2 = getImageID();
	roiManager("select",0);		
	run("Enlarge...", "enlarge=-20");
	run("Enlarge...", "enlarge=20");
	setColor(0);
	fill();
	setThreshold(0, 0);	
	run("Analyze Particles...", "size=10-Infinity add");
	roiManager("select", 0);
	roiManager("delete");
	keepLargestROI();
	roiManager("select", 0);	
	roiManager("update");
	run("Select None");
	selectImage(id2); close();	
	run("Select None");
}

function keepLargestROI() {
	/* keep the largest ROI from the ROI manager */
	n = roiManager("count");
	run("Select None");
	amax = 0;
	istar = -1;
	for (i = 0; i < n; i++) {
		roiManager("select",i);
		a = getValue("Area");
		if (a > amax) {
			istar = i;
			amax = a;
		}
	}
	for (i = n-1; i >= 0; i--) {		
		if (i != istar) {
			roiManager("select",i);
			roiManager("delete");
		}
	}	
	roiManager("select",0);	
}

function EMThreshold() {
	/* 
	 *  Fit a Gaussian mixture with 2 components using an 
	 * 	expectation minimization algorithm to define a threshold
	 * 	
	 */	
	id = getImageID();	
	N = getWidth() * getHeight();
	getHistogram(values, counts, 255);	
	n = values.length;
	
	getStatistics(area, mean, min, max, std, histogram);
	m1 = (mean-4*std);
	s1 = std*std/4;
	a1 = 0.5;
	
	m2 = (mean+4*std);
	s2 = std*std/4;
	a2 = 0.5;
	
	p1 = newArray(n);
	p2 = newArray(n);	
	y = newArray(n);
	for (iter = 0; iter < 100; iter++) {			
		// E
		for (i = 0; i < n; i++) {			
			y1 = gaussian1(values[i],a1,m1,s1);
			y2 = gaussian1(values[i],a2,m2,s2);
			y = y1 + y2;
			p1[i] = y1 / (y1 + y2);
			p2[i] = y2 / (y1 + y2);			
		}
		// M	
		m1=0;s1=0;a1=0;
		m2=0;s2=0;a2=0;
		N = 0;
		for (i = 0; i < n; i++) {	
			m1 += p1[i] * values[i] * counts[i];
			s1 += p1[i] * values[i] * values[i] * counts[i];
			a1 += p1[i] * counts[i];			
			
			m2 += p2[i] * values[i] * counts[i];
			s2 += p2[i] * values[i] * values[i] * counts[i];
			a2 += p2[i] * counts[i];
			N += counts[i];
		}
		
		m1 /= a1;
		s1 /= a1; s1 -= m1*m1;

		m2 /= a2;
		s2 /= a2; s2 -= m2*m2;
				
		a1 = a1 / N;
		a2 = a2 / N;		
	}		
	// likelihood test
	t = mean + 3*std;
	for (i = 0; i < values.length-1; i++) {
		if (p1[i] > 0.5 && p1[i+1] <= 0.5){
			t = values[i];
		}
	}
	setThreshold(t,max);
}

function gaussian1(x,a,m,s) {	
	d = abs(x-m);
	p = a * exp(-0.5*d*d/s) / sqrt(2*PI*s);	
	return p;
}

function gradTheCurve(x,y,tmax,Imax) {
	// find the contour
	dx = diff(x);
	dy = diff(y);
	delta = 0; // displacement of the curve
	K = 4*tmax;	
	for (i = 0; i < x.length; i++) {
		n = sqrt(dx[i]*dx[i]+dy[i]*dy[i]);
		swx = 0;
		sw = 0;		
		m = 0;
		maxval = 0;
		xmax = x[i];
		ymax = y[i];
		for (k = 0; k < K; k++) {
			t = -tmax + 2 * tmax * k / K;
			xt = x[i] + t * dy[i] / n;
			yt = y[i] - t * dx[i] / n;
			I = getPixel(xt,yt);
			if (I>maxval) {
				maxval = I;
				xmax = xt;
				ymax = yt;
			}
			w = I * exp(-0.5*(t*t)/(tmax*tmax));
			swx += w * t;
			sw += w;
			m++;
		}
		//x[i] = xmax;
		//y[i] = ymax;
		if (sw > 0) {			
			tstar = swx / sw;
			vx = tstar * dy[i] / n;
			vy = - tstar * dx[i] / n;
			V = sqrt(vx*vx+vy*vy);	
			r = 1 / ( 1 + exp(-20*sw/Imax/m));			
			x[i] = x[i] + r*vx;
			y[i] = y[i] + r*vy;					
			delta = delta +  V;
		}
		Imax = maxval;
	}	
	return delta / x.length;
}

function measureIntensityOutside(id) {
	/* Measure intensity outside the organoid
	 *  
	 * Use a simple threshold to discard the darker organoid and shrink the selection
	 * to avoid taking into account regions from inside the organoid. return the mean intensity
	 */
	selectImage(id);	
	run("Select None");
	setAutoThreshold("Minimum dark");
	run("Create Selection");	
	run("Enlarge...", "enlarge=-50");
	val = getValue("Mean");
	run("Select None");	
	resetThreshold;
	return val;
}

function transparency() {
	// Measure the transparency as the Mean of a filtered image
	run("Duplicate...", " ");
	run("32-bit");
	run("Gaussian Blur...", "sigma=1");
	run("Convolve...", "text1=[-1 -1 -1\n-1 8 -1\n-1 -1 -1\n] normalize");
	run("Abs");
	v = getValue("Mean");
	close();
	return v;
}

function getContourInflectionPoints(curvature,r0) {
	// return list of index where the curvature is more than 2 time the minor axis.
	//r0 = getValue("Major");
	tmp = smooth1d(curvature,2);
	//tmp = curvature;	
	threshold = 1 / r0;
	inflx = findPeaks(tmp, threshold, 3);
	
	xp = newArray(inflx.length);
	yp = newArray(inflx.length);
	for (i=0;i<inflx.length;i++){
		xp[i] = inflx[i];
		yp[i] = tmp[inflx[i]]; 
	}	
	
	id = getImageID();
	/*
	Plot.create("Title", "X","Y");
	Plot.add("line", Array.getSequence(curvature.length), curvature);
	Plot.setColor("green");
	Plot.add("line", Array.getSequence(curvature.length), tmp);	
	xl = newArray(0,curvature.length-1);
	yl = newArray(threshold,threshold);
	Plot.setColor("red");
	Plot.add("line",xl,yl);
	Plot.add("circle",xp,yp);
	Plot.update();
	*/
	selectImage(id);
	
	return inflx;
}

function findPeaks(A,t,p) {
	/* find local maxima in neighborhood [-p,p] above threshold */
	maximas = newArray(A.length);
	k = 0;
	for (i = 0; i < A.length; i++) if (A[i] > t) {			
		flag = true;
		for (j = i-p; j <= i+p && flag; j++) {			
			if (A[getCircBB(A.length,j)] > A[i]) {
				flag = false;				
			}
		}
		if (flag) {
			maximas[k] = i;
			k++;
		}
	}	
	return Array.trim(maximas, k);
}

function getCircBB(n,i) {
	/* index with circular boundary conditions */
	if (i >= 0 && i < n) return i;
	if (i < 0) return n+i;
	if (i >=  n) return i-n;
}

function getCurvilinearLength(x,y,pixel_size) {
	// return the distance between consecutive points using circular boundary conditions	
	s = newArray(x.length);	
	for (i = 0; i < x.length-1; i++) {
		dx = x[i+1] - x[i];
		dy = y[i+1] - y[i];
		s[i] = pixel_size * sqrt(dx*dx+dy*dy);
	}
	dx = x[0] - x[x.length-1];
	dy = y[0] - y[y.length-1];
	s[x.length-1] = pixel_size * sqrt(dx*dx+dy*dy);
	return s;
}

function getWeightedEnergy(value,weight) {
	// compute the energy of the curvature (sum weigth * value^2) / sum(weight)
	e = 0;
	p = 0;
	for (i = 0; i < value.length; i++) {
		e = e + value[i] * value[i] * weight[i];
		p = p + weight[i];
	}
	if (p>0) {return e / p;}
	else {return 0;}
}

function DNDhelper(x0,y0,x1,y1,x2,y2,pixel_size) {
	// Norm of the variation of the normal projected on the tangent times 
	// compute derivatives at point p-1 and p+1
	dx1 = (x1 - x0) * pixel_size;
	dx2 = (x2 - x1) * pixel_size;
	dy1 = (y1 - y0) * pixel_size;
	dy2 = (y2 - y1) * pixel_size;
	s1 = sqrt(dx1 * dx1 + dy1 * dy1);
	s2 = sqrt(dx2 * dx2 + dy2 * dy2);
	dx1 = dx1 / s1;
	dx2 = dx2 / s2;
	dy1 = dy1 / s1;
	dy2 = dy2 / s2;
	// derivative of the normal n=(dy,-dx) 
	dnx = (dy2-dy1)/s1;
	dny = -(dx2-dx1)/s1;
	// project on tangent u=(dx,dy): dn.u
	dnu = dnx * dx1 + dny * dy2;
	return sqrt(dnu*dnu);
}

function DND(x,y,pixel_size) {
	n = x.length;
	e = newArray(n);
	for (i = 1; i < n-1; i++) {
		e[i] = DNDhelper(x[i-1],y[i-1],x[i],y[i],x[i+1],y[i+1],pixel_size);
	}
	e[0] = DNDhelper(x[n-1],y[n-1],x[0],y[0],x[1],y[1],pixel_size);
	e[n-1] = DNDhelper(x[n-2],y[n-2],x[n-1],y[n-1],x[0],y[0],pixel_size);
	return e;
}

/* Display of contours with color coding */

function colorContour(x,y,value,width) {
	// color the contour with value as overlay
	lut = getLutHexCodes("Ice");		
	if (colorbarmin==0 && colorbarmax==0) {
		Array.getStatistics(value, colorbarmin, colorbarmax, mean, stdDev);
		colorbarmin = - maxOf(abs(colorbarmax), abs(colorbarmin));
		colorbarmax = -colorbarmin;
	}
	for (i = 0; i < x.length-1; i++) {
		cidx = Math.floor((lut.length-1)*(value[i] - colorbarmin) / (colorbarmax-colorbarmin));
		cidx = minOf(lut.length-1, maxOf(0, cidx));
		setColor(lut[cidx]);
		Overlay.drawLine(x[i], y[i], x[i+1], y[i+1]);
		Overlay.show;
		Overlay.setStrokeWidth(width);		
	}
	// close the  contour
	cidx = Math.floor((lut.length-1)*(value[i] - colorbarmin) / (colorbarmax-colorbarmin));
	cidx = minOf(lut.length-1, maxOf(0, cidx));	
	setColor(lut[cidx]);
	Overlay.drawLine(x[x.length-1], y[y.length-1], x[0], y[0]);
	Overlay.show;
	run("Labels...", "color=white font=12");
}

function getLutHexCodes(name) {
	//setBatchMode("hide");
	newImage("stamp", "8-bit White", 10, 10, 1);
   	run(name);
   	getLut(r, g, b);
   	close();
	dst = newArray(r.length);
	for (i = 0; i < r.length; i++) {
		dst[i] = "#"+dec2hex(r[i]) + dec2hex(g[i]) + dec2hex(b[i]);
	}   
	//setBatchMode("show");
	return dst;	
}

function dec2hex(n) {
	// convert a decimal number to its heximal representation
	codes = newArray("0","1","2","3","4","5","6","7","8","9","A","B","C","D","E","F");
	k1 = Math.floor(n/16);
	k2 = n - 16 * Math.floor(n/16);
	return codes[k1]+codes[k2]; 
}

/* Contour curvature */

function getAverageRadius(x,y,pixel_size) {
	/* compute the average radius of the polygone defined by x,y */
	n = x.length;
	xbar = 0;
	ybar = 0;
	for (i = 0; i < n; i++) {
		xbar += x[i];
		ybar += y[i];
	}
	xbar /= n;
	ybar /= n;
	r = 0;
	for (i = 0; i < n; i++) {
		dx = x[i]-xbar;
		dy = y[i]-ybar;
		r += sqrt(dx*dx+dy*dy);
	}		
	return r/n * pixel_size;
}

function getContourCurvatureGeo(x,y,pixel_size) {
	// Compute cuvature with periodic boundary conditions using a geometric method
	n = x.length;
	kappa = newArray(n);
	for (i = 0; i < n; i++) {
		ip = getCircBB(n,i-1);
		in = getCircBB(n,i+1);		
		kappa[i] = curvHelper(x[ip],y[ip],x[i],y[i],x[in],y[in]) / pixel_size;		
	}
	return kappa;
}

function curvHelper(x0,y0,x1,y1,x2,y2) {
	/* Use Heron formula to compute curvature
	 as the radius of the osculating circle
	 using a geometric approach
	 https://scholar.rose-hulman.edu/cgi/viewcontent.cgi?article=1233&context=rhumj
	*/
	xa = x0 - x1;
	ya = y0 - y1;
	
	xb = x2 - x1;
	yb = y2 - y1;
	
	xc = x2 - x0;
	yc = y2 - y0;
	
	a = sqrt(xa*xa+ya*ya);
	b = sqrt(xb*xb+yb*yb);
	c = sqrt(xc*xc+yc*yc);
	//S = sqrt((a+b+c)*(-a+b+c)*(a-b+c)*(a+b-c)) / 4;
	//R = a*b*c/(2*S);	
	Delta = (xb*ya-xa*yb) / 2.0;
	Kappa = 4.0 * Delta / (a*b*c);
	return Kappa;	
}


function getContourCurvature(x,y) {
	/* Compute cuvature with periodic boundary conditions */
	n = x.length;
	gamma = newArray(n);
	for (i = 1; i < n-1; i++) {
		dx = (x[i+1] - x[i-1])/2;
		dy = (y[i+1] - y[i-1])/2;
		dxx = x[i-1] - 2 * x[i] + x[i+1];
		dyy = y[i-1] - 2 * y[i] + y[i+1];
		gamma[i] = (dx * dyy - dy * dxx) / pow(dx*dx + dy*dy, 3.0/2.0);
	}
	// at 0 with perdiodic bc
	dx = x[1] - x[n-1];
	dy = y[1] - y[n-1];
	dxx = - x[n-1] + 2 * x[0] - x[1];
	dyy = - y[n-1] + 2 * y[0] - y[1];
	gamma[0] = (dx * dyy - dy * dxx) / pow(dx*dx + dy*dy, 3.0/2.0);
	// at n-1 with periodic bc
	dx = x[n-2] - x[0];
	dy = y[n-2] - y[0];
	dxx = - x[n-2] + 2 * x[n-1] - x[0];
	dyy = - y[n-2] + 2 * y[n-1] - y[0];
	gamma[n-1] = (dx * dyy - dy * dxx) / pow(dx*dx + dy*dy, 3.0/2.0);
	return gamma;
}

function diff(u) {
	/* First derivatives with periodic bounday conditions */
	n = u.length;
	d = newArray(n);
	for (i = 0; i < n - 1; i++) {
		d[i] = u[i+1] - u[i];
	}
	d[n-1] = u[0] - u[n-1];
	return d; 
}

function diff2(u) {
	/* Second derivatives with periodic bounday conditions */
	n = u.length;
	d = newArray(n);
	for (i = 1; i < n-1; i++) {
		d[i] = - u[i-1] + 2 * u[i] - u[i+1];
	}
	d[n-1] = -u[n-2] + 2 * u[n-1] - u[0];
	d[0]   = -u[n-1] + 2 * u[0]   - u[1];  
	return d; 
}

/* FFT */

function roiFractalDimension(x, y, doplot) {
	/* Compute the fractal dimension of the contour
	 * should be between 1 and 2..
	 * https://arxiv.org/abs/1201.3097		
	 */
	 		
	if (x.length < 8) return 1.0;
			
	if (doplot == true) {
		Plot.create("curve","x","y",x,y);
	}
	
	m = powerspectrum(x, y);
	
	f = newArray(m.length-1);
	mf = newArray(m.length-1);
	for (i = 0; i <f.length; i++) {	
		f[i] = log(i+1);
		mf[i] = log(m[i+1]);
	}
	
	if (doplot == true) {
		Plot.create("power spectrum","k","psd",f,mf);
	}
	
	Fit.doFit("Straight Line", f, mf);
	if (doplot == true) {
		Fit.plot;
	}
	
	beta = Fit.p(1);
	return 1-(beta + 2) / 6;	
}

function powerspectrum(_real,_imag) {
	/* Log power spectrum*/
	N = nextPowerOfTwo(_real.length);
	//print(N);
	real = Array.resample(_real,N);
	imag = Array.resample(_imag,N);
	//real = padarray(_real,N);
	//imag = padarray(_imag,N);
	fft(real,imag);
	m = newArray(real.length/2);
	for (i = 0; i < m.length; i++) {
		m[i] = real[i]*real[i]+imag[i]*imag[i];
	}
	return m;
}

function roiFftSmooth(_threshold) {
	// Smooth the contour using thresholding of the FFT
	run("Interpolate", "interval=5 smooth adjust");
	//Roi.getCoordinates(x, y);	
	Roi.getSplineAnchors(x, y);	
	fftSmooth(x,y,_threshold);
	setColor("green");
	//makeSelection("freehand", x, y);
	Roi.setPolygonSplineAnchors(x, y);
	run("Interpolate", "interval=10 smooth adjust");
}

function fftSmooth(_real,_imag,_threshold) {
	//  Smoothing by Fourier coefficient hard thresholding
	N = nextPowerOfTwo(_real.length);	
	real = padarray(_real,N);
	imag = padarray(_imag,N);
	fft(real,imag);	
	for (i = 2; i < real.length; i++) {
		m = sqrt(real[i]*real[i]+imag[i]*imag[i]);	
		if (m < _threshold*_real.length) {
			real[i] = 0;
			imag[i] = 0;
		} 
	}	
	ifft(real,imag);	
	for (i = 0; i < _real.length; i++) {
		_real[i] = real[i];
		_imag[i] = imag[i];
	}		
}

function nextPowerOfTwo(_N) {
	// Return the next power of two integer
	return pow(2.0, 1+floor(log(_N)/log(2)));
}

function padarray(_src,_N) {
	//  Extend the size of the array to size N and fill it with zeros 
	if (_N > _src.length) {
	dest = newArray(_N);
	for (i = 0; i < _src.length; i++) {
		dest[i] = _src[i];
	}
	} else {
		dest=_src;
	}
	return dest;
}

function fft(_real,_imag) {
	// Forward Fast fourier transform
	_fft_helper(_real,_imag,0,_real.length,true);
}

function ifft(_real,_imag) {
	// Inverse fast fourier transform
	_fft_helper(_real,_imag,0,_real.length,false);
	for(k = 0; k < _real.length; k++) {
		_real[k] = _real[k] / _real.length;
		_imag[k] = _imag[k] / _imag.length;
	}
}

function _fft_helper(ax,ay,p,N,forward) {
	// Fast Fourier Transform helper function
    if(N < 2) {
    } else {
        separate(ax,p,N);
        separate(ay,p,N);        
        _fft_helper(ax, ay, p,       N/2, forward);  
        _fft_helper(ax, ay, p + N/2, N/2, forward);                               
        for(k = 0; k < N/2; k++) {
        	// even
            ex = ax[k    +p];
            ey = ay[k    +p];
            // odd
            ox = ax[k+N/2+p];
            oy = ay[k+N/2+p];
            // twiddle factor
            if (forward == true) {
            	wx = cos(-2*PI*k/N);
            	wy = sin(-2*PI*k/N);
            } else {
            	wx = cos(2*PI*k/N);
            	wy = sin(2*PI*k/N);
            }
            sx = wx * ox - wy * oy;
            sy = wx * oy + wy * ox;
            ax[k    +p] = ex + sx;
            ay[k    +p] = ey + sy;
            ax[k+N/2+p] = ex - sx;
            ay[k+N/2+p] = ey - sy;                      
        }           
    }
}

function separate(a, p, n) {
	// 	Separate odd and even
	b = newArray(n/2);	
    for(i=0; i < n/2; i++) {
        b[i] = a[2*i+1+p];        
    }
    for(i=0; i < n/2; i++) {
        a[i+p] = a[2*i+p];        
    }
    for(i=0; i < n/2; i++) {
        a[i+n/2+p] = b[i];        
    }
}

/*** END oF FFT **/

function smooth1d(src,sigma) {
	n = src.length;
	dst = newArray(n);	
	k = Math.ceil(3*sigma);
	for (i = 0; i < n; i++) {
		sum = 0;
		count = 0;
		for (j = i-k; j <= i+k; j++) {
			d = (i-j)/sigma;
			w = exp(-0.5*d*d);
			sum =  sum + w * src[getCircBB(n,j)];
			count = count + w;			
		}
		dst[i] = sum / count;
	}
	return dst;
}

function closeRoiManager() {
	/* Empty the ROI manager and close it*/
	while (roiManager("count")>0) {roiManager("select", roiManager("count")-1);roiManager("delete");}
	close("ROI Manager");
	run("Select None");
}

function main() {
	/* entry point */
	if (!endsWith(filename,".tif")) {
		exit();
	}
	start_time = getTime();
	run("Close All");
	run("Set Measurements...", "  redirect=None decimal=10");
	closeRoiManager();
	print("-------------------------------------------------");
	print(filename);
	folder = File.getDirectory(filename);
	fname = File.getName(filename);
	tbl = "Organoid Morphology";
	n = initTable(tbl);
	processFile(n, folder, fname, tbl);	
	run("Collect Garbage");	
	end_time = getTime();
	print("Elapsed time " + (end_time - start_time) / 1000 + " seconds.");
}

main();
