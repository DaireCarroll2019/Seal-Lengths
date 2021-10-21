# Seal-Lengths
A custom function to measure the lengths and widths of polygons representing seals automatically identified from drone images

Seal_Length_Function contains the function to extract the curved length and width from a sf spatial polygon (assuming the correct projection system, values will be given in the approprioate projection unit)

Seal_Length_Process contains a pipeline to extract and annotate seal lengths and widths from a folder containing multiply .shp file subfolders. Each .shp file subfolder will be opened, polygons measured, and the .shp file overwritten with a new file containing infomation about polygon dimensions based on the Seal_Length_Function. The pipeline contains an interactive step during which users are required to view polygons above 1.6 m in length and determine wether they should remain in the file or be removed.
