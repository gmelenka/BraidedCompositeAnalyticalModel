# BraidedCompositeAnalyticalModel
![Braided Composite Analytical Model](https://raw.githubusercontent.com/gmelenka/BraidedCompositeAnalyticalModel/master/BraidModelLogo-01_256.png)

Analytical Model for analyzing the mechanical properties of tubular braided composites

Braided Composites can be manufactured in several different braiding patterns: Diamond (1/1), Regular (2/2) and Hercules (3/3).  Examples of the different patterns are shown below

![Diamond (1/1)](https://raw.githubusercontent.com/gmelenka/BraidedCompositeAnalyticalModel/master/diamond.jpeg)

![Regular (2/2)](https://raw.githubusercontent.com/gmelenka/BraidedCompositeAnalyticalModel/master/regular.jpeg)

![Hercules (3/3](https://raw.githubusercontent.com/gmelenka/BraidedCompositeAnalyticalModel/master/hercules.jpeg)

Current work on braided composites typically focuses on one specific pattern.  The goal of this model is to be able to predict the mechanical properties of all the available braiding patterns.
Another limitation of current studies that examine braided composites is that most researchers do not take into consideration the geometry and the physical process of manufacturing composite braids.
Both braid geometry and the manufacturing process has an effect on the mechanical properties of tubular braided composites.

![Braided Composite AAnalytical Model ScreenShot](https://raw.githubusercontent.com/gmelenka/BraidedCompositeAnalyticalModel/master/braidedCompositeModulus.png)


Installation:

The program can be installed using several approaches

1) Install as a Matlab App
The Braided Composite Analytical Model has been developed using the Matlab App Developer program which is available in Matlab 2016a
The Braided Composite Analytical Model can be installed by selecting the file Braided Composite Elastic Properties.mlappinstall
This will add the Braided Composite Analytical Model to the Apps ribbon bar in Matlab

2) Install as Exe
The Braided Composite Analytical Model can be run as an executable file.
Select the file BraidedCompositeAnalyticalModel.exe to install the program. A MATLAB runtime engine will be required to run this file.  The runtime engine will automatically download the first timee the executable is run.
Once the runtime engine has been isntalled the file BraidedCompositeAnalyticalModel.exe can be used.  Once installed, this option will not require MATLAB.

3) Run MATLAB scripts
The main braid model has been written as matlab m-file functions.  The m-file functions can be used within the MATLAB environment without the need of a graphical user interface.  In order to run the braid model the
function braidModel.m will be used.  The input and output of the braidModel function is shown below.  The function will return the mechanical properties of a composite braid.

[Ex, Ey, Ez, GxyCombined, GyzCombined, GzxCombined] = braidModel(S, Sm, angle, n, r0, a, b, beta, braidType)