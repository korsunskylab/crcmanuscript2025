Each file contains an SFC object where
- column 1: serial number of FOV
- column 2: sfc POLYGON containing the bounding box for that FOV. This is calculated from the detected transcripts provided by Vizgen. Each tx was assigned to an FOV, so I can calculate the bounding box for an FOV by drawing a bbox around all tx assigned to that FOV.
-> convenience for plotting

  