This package includes Matlab scripts to analyze the morphology of volcanic edifices using two main algorithms - DrianageVolc and MorVolc. DrianageVolc is designed to analyze drainage basins, networks, and divides of volcanoes. MorVolc is a recreation of the original IDL script originally published by Grosse et al. (2012) to analyze volcanic edifice geometries. Complete descriptions of each algorithm can be found in the DrainageVolc_Analysis.m and MorVolc_Analysis.m files, respectively. Furthermore, example input scripts can be found in the Example_Scripts folder.

The paper that introduces these algorithms is currently submitted to Geology, and should be referenced as:

O'Hara, D., Goren, L., van Wees, R.M.J., Campforts, B., Grosse ,P., Lahitte, P., Kereszturi, G., Kervyn, M. (submitted). Volcano drainage morphology co-varies with age of activity - new insights on radial drainage development and edifice erosion. Geology.

Both algorithms use TopoToolbox for their analysis, which can be downloaded from https://topotoolbox.wordpress.com/ and placed into the General_Scripts folder. Other community-written scripts are also used by the algorithms to perform simple tasks. These scripts are already incorporated into General_Scripts folder, following CC licencing. Script references are provided below.

fitellipse:
Fitzgibbon, Andrew W., Maurizio Pilu, and Robert B. Fisher. "Direct least squares fitting of ellipses". IEEE Transactions on Pattern Analysis and Machine Intelligence, 21(5), 476-480, May 1999. http://homepages.inf.ed.ac.uk/rbf/CVonline/LOCAL_COPIES/FITZGIBBON/ELLIPSE/

idw:
Andres Tovar (2022). Inverse distance weight function (https://www.mathworks.com/matlabcentral/fileexchange/46350-inverse-distance-weight-function), MATLAB Central File Exchange. Retrieved December 28, 2022

ll2utm:
Francois Beauducel (2015). https://github.com/IPGP/mapping-lib/blob/master/latlonutm/ll2utm.m

viridis:
Ander Biguri (2022). Perceptually uniform colormaps (https://www.mathworks.com/matlabcentral/fileexchange/51986-perceptually-uniform-colormaps), MATLAB Central File Exchange. Retrieved December 28, 2022.

Kringing:
Wolfgang Schwangart (2010). https://github.com/wschwanghart/kriging

Other Reference:
Grosse, P., De Vries, B. V. W., Euillades, P. A., Kervyn, M., & Petrinovic, I. A. (2012). Systematic morphometric characterization of volcanic edifices using digital elevation models. Geomorphology, 136(1), 114-131.
