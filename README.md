#SegmentationCheck
This code is developed to analyze XMT (X-ray microtomography) images from aluminium copper samples. Fig1 shows a reconstruted volume from a set of XMT images. 
<p align="center"><img src=images/img1.png width="500"></p>
<p align="center"><i>fig1.</i> Reconstruced volume from XMT images.</p>
##Code Description
First section of the code, `comp()`, is used  to convert mass fraction to volume fraction in eutectic alloys, using either Schiel equation or lever rule depending on user's intentions.
Second section, `VolFracTomog()`, is used to calculate the area fraction of the secondary phase in a micrograph. User can implement K-means method, use their own threshold array or use Ostu's method to segment the image. The code then can be used to iterate over a number of images, this is specially useful if the user is inteding to determine the volume fraction of a phase in XMT images. 
The thir bit of the code is a custom code to generate a plot for a specific problem, user may write their own code to visualize their results.
##Sample results
For this case, I have run this code for 4 different materials. Note that close to 11000 XMT images (fig2) were analyzed in the process. 
<p align="center"><img src=images/img2.png width="500"></p>
<p align="center"><i>fig2.</i> Sample XMT image.</p>
A sample result image has also been generated (fig3). Note that this is only a sample illustration and does not represent the actual results generated from the XMT images.
<p align="center"><img src=images/img3.png width="500"></p>
<p align="center"><i>fig3.</i> Sample illustartion of the results.</p>
# Tested versions
These scripts have been tested using:
- Python 2.7.11 Anaconda 2.3.0
