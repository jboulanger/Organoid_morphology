# Organoid Morphology

This repository contains an ImageJ 1 macro that analyses the morphology of brain organoids slices acquired in brightfiled microscopy.

![B8-kit-day18-2](https://github.com/jboulanger/Organoid_morphology/assets/3415561/1be0f509-007b-4457-966d-0f7a09be3281)

The segmented contour is analyzed to extract several parameters: 
- Area
- Perimeter
- Average Radius
- Roundess
- Aspect Ratio
- Feret
- min Feret
- Circularity
- number of inflection points
- sum of weighted curvature
- DNE
- estimated transparency
- mean of the curvature
- standard deviation of the curvature
- Average radius time the standard deviation of the curvature.

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


