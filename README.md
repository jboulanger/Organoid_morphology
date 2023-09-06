# Organoid Morphology

This repository contains an ImageJ 1 macro that analyses the morphology of brain organoids slices acquired in brightfiled microscopy.

![B8-kit-day18-2](https://github.com/jboulanger/Organoid_morphology/assets/3415561/1be0f509-007b-4457-966d-0f7a09be3281)

The segmented contour is analyzed to extract several parameters: 
- Area [um^2]: Area of the selection in um^2 ("Area" measurement)
- Perimeter [um]: Perimeter of the selection ("Perim." measurement)
- Average Radius [um]: Average distance R0 of the contour to the center of the selection
- Roundess: 4 x [Area/π (Major axis)2] ("Round" measurement) 
- Aspect Ratio: Major axis / Minor axis of the fitted ellipse ("AR" measurement)
- Feret [um]: the longest distance between any two points along the selection boundary ("Feret" Measurement)
- Min Feret [um]: the minimum distance between any two points along the selection boundary ("MinFeret" Measurement)
- Circularity: 4π x (Area/Perimeter^2) ("Circ." Measurement)
- Inflection points: number of detected inflection points 
- Weighted curvature [um^-1] : Logarithm of the sum of the curvature x segment length
- DNE: logarithm of the square of the variation of the normal n=(dy,-dx) of the contour projected on its tangent t=(dx,dy) where dx and dy are the first derivative in x and y. DNE is normalised by the average radius (R0).
- Transparency: the mean response of the Laplacian of Gaussian (LoG) filter.
- Mean curvature [um^-1]: Average of the curvature along the contour
- Std curvature: Standard deviation of the curvature along the contour
- R0 x Std curvature: Standard deviation of the curvature along the contour normalized by the average radius (R0)

The curvature along the contour is the inverse of the radius of the osculating circle and is computed as [(dx * dyy) – (dy * dxx)] / [(dx^2) + (dy^2)]^3/2 where dx and dy are the first derivative in x and y and dxx and dyy are the second derivative of the contour. The curvature is computed with a [geometric approach](https://scholar.rose-hulman.edu/cgi/viewcontent.cgi?article=1233&context=rhumj) using Heron's formula.



## Installation
Download the macro [Organoid_Morphology.ijm](https://raw.githubusercontent.com/jboulanger/Organoid_morphology/main/Organoid_Morphology.ijm).

## Usage
Open the macro in the script editor and press either Run or Batch. 

- The input file must be a TIFF file with a valid pixel calibration in microns.
- Select manual selection if a manual segmentation is preferred
- Overlay: select the measure to overlay with the image
- Display Info: select this to display additional information as an overlay on the image
- Add colorbar: select this to display a color bar to the selected measure
- Colorbar min: set the minimum value of the colorbar
- Colorbar max: set the maximum value of the colorbar
- Save image as jpeg and close: select this to save the image as a jpg file with the annotations (contour, info, colorbar, etc)
- Use Saved ROI: select this to use previoulsy stored ROI (filename.zip -> filename.roi)
- Marker scale: define the scale of the marker for the inflection points

Measurements are appended into a single table.

A typical workflow would require to process many images running the macro using "Batch" and saving jpg and roi files.
Inspect the jpeg file to identify the poorly segmented file and correct those with a manual segmentation step to update the jpeg and roi file. 
Finally, run again the macro on all files with the option "Use Saved ROI" enabled.


