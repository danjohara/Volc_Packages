This package includes Matlab scripts to analyze the morphology of volcanic edifices using two main algorithms - DrianageVolc and MorVolc. DrianageVolc is designed to analyze drainage basins, networks, and divides of volcanoes. MorVolc is a recreation of the original IDL script originally published by Grosse et al. (2012) to analyze volcanic edifice geometries. Complete descriptions of each algorithm can be found in the DrainageVolc_Analysis.m and MorVolc_Analysis.m files, respectively. Furthermore, example input scripts can be found in the Example_Scripts folder.

The paper that introduces these algorithms is in process to be submitted to Earth Surface Dynamics, and should be referenced as:

O'Hara, D., Goren, L., van Wees, R.M.J., Campforts, B., Grosse ,P., Lahitte, P., Kereszturi, G., Kervyn, M. (submitted). Volcano drainage morphology co-varies with age of activity - new insights on radial drainage development and edifice erosion. Earth Surface Dynamics.

Both algorithms use TopoToolbox for their analysis, which can be downloaded from https://topotoolbox.wordpress.com/ and placed into the General_Scripts folder. These algorithms also use Matlab's built-in packages that need downloaded; including Statistics and Machine Learning, Parallel Computing, Image Processing, Mapping, and Signal Processing. Other community-written scripts are also used by the algorithms to perform simple tasks. These scripts are already incorporated into General_Scripts folder, following CC licencing. Script references are provided below.

fitellipse:
Fitzgibbon, Andrew W., Maurizio Pilu, and Robert B. Fisher. "Direct least squares fitting of ellipses". IEEE Transactions on Pattern Analysis and Machine Intelligence, 21(5), 476-480, May 1999. http://homepages.inf.ed.ac.uk/rbf/CVonline/LOCAL_COPIES/FITZGIBBON/ELLIPSE/

idw:
Andres Tovar (2022). Inverse distance weight function (https://www.mathworks.com/matlabcentral/fileexchange/46350-inverse-distance-weight-function), MATLAB Central File Exchange. Retrieved December 28, 2022

ll2utm:
Francois Beauducel (2015). https://github.com/IPGP/mapping-lib/blob/master/latlonutm/ll2utm.m

viridis:
Ander Biguri (2022). Perceptually uniform colormaps (https://www.mathworks.com/matlabcentral/fileexchange/51986-perceptually-uniform-colormaps), MATLAB Central File Exchange. Retrieved December 28, 2022.

bluewhitered:
Nathan Childress (2023). bluewhitered (https://www.mathworks.com/matlabcentral/fileexchange/4058-bluewhitered), MATLAB Central File Exchange. Retrieved May 15, 2023.

Kringing:
Wolfgang Schwangart (2010). https://github.com/wschwanghart/kriging

