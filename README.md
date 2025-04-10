1. Operating System:
Windows 10/11, MATLAB

2. File Format Requirements:
 The input file must be a .txt file. Each column of data should be separated by commas.
 The first column: X-coordinate of the observation point
 The second column: Y-coordinate of the observation point
 The third column: Observation height
 The fourth column: Observed value
 The fifth column: Weight

4. Code Description:
(1) "Focusing_Inversion" is the main program. Simply run this file to start. The interface includes parameter settings and visualization of observation data.The main GUI interface of the software consists of three sections. Section one is the “Parameters Settings,” which includes the focusing factor, the iteration parameters such as maximum iterations, the adaptive factor and tolerance, and the depth weighting parameters. Additionally, there is a checkbox for selecting constraints, which includes depth weighting, data weighting, reference model, and solution bounds. Section two is the “Subsurface Model,” which defines the parameters for the x, y, and z axis ranges. If the reference model checkbox is checked, additional parameters for the reference model are also provided. Section three is the “Observation Data.” In this section, we can select which gravity component to use as the observation data. A button is provided to load the corresponding gravity grid data. After data loading is complete, there is another button to display the map of the input data.
(2) "Result" is the result display interface, which runs automatically after the computation is complete.the results GUI interface of the software consists of two parts. Part one displays a table of the 3D density solution along with statistics of the solution. Additionally, there is a subsection for result visualization, where the solution slice can be shown along the x, y, and z axes at the defined positions. Part two provides functionality for using the 3D density solution to predict gravity components, compute the forward modeling results’ statistics, and display the corresponding maps.
