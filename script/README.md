Example scripts for simulation studies and data analysis used in the paper.

1. segmentation_test.py

   Segmentation simulation study.

2. geopip_test.py

   GeoPIP simulation study.

   Simulate data based on the GeoPIP model, and make inference using PhyML, CTMC+NJ, PIP+NJ, GeoPIP+NJ.

3. hpip_test.py

   hPIP simulation study.

   Simulate date based on the hPIP model, and make inference using PhyML, CTMC+NJ, PIP+NJ, GeoPIP+NJ.

4. indelible_test.py

   INDELible simulation study.

   Simulate date using INDELible, and make inference using PhyML, CTMC+NJ, PIP+NJ, GeoPIP+NJ.

   Note: this file should be used in pair with “control.txt”, which sets parameters for running INDELible. Different simulation setups can be used by changing BOTH the “control.txt” file and relevant parameters in the first part of the “indelible_test.py” file.

5. muscle_test.py

   MUSCLE simulation study.

   Simulate data using INDELible, estimate alignment using MUSCLE, and make inference on the aligned data using PhyML, CTMC+NJ, PIP+NJ, GeoPIP+NJ.

   Note: this file should be used in pair with “control.txt”, which sets parameters for running INDELible. Different simulation setups can be used by changing BOTH the “control.txt” file and relevant parameters in the first part of the “muscle_test.py” file.

6. data_analysis.py

   Data analysis using the Molluscan data in the ../data folder.

   Data analysis using PhyML, CTMC+NJ, PIP+NJ, GeoPIP+NJ.

   Note: to analyze a new data set, just change the input file location in “data_analysis.py” to the location of the new fasta file.
