import csv, sys
import os

# Python Script to Convert GSLIB Geo-EAS to R comma delimited
# python GSLIB2R.py input_filename output_filename
# Author Michael Pyrcz, University of Texas at Austin, @GeostatsGuy

# Convert GSLIB Geo-EAS files to a pandas DataFrame for use with Python methods
def GSLIB2Dataframe(data_file):
    import os
    from pathlib import Path
    import numpy as np  
    import pandas as pd

    colArray = []
# Open file, exit if file doesn't exit    
    my_file = Path(data_file)
    if my_file.is_file() == False:
    	print("Error - file not found!")
    	quit()
    	
# Read file GSLIB file into DataFrame    
    with open(data_file) as myfile:   # read first two lines
        head = [next(myfile) for x in range(2)]
        line2 = head[1].split()
        line2
        ncol = int(line2[0])
        for icol in range(0, ncol):
            head = [next(myfile) for x in range(1)]
            colArray.append(head[0].split()[0])
        data = np.loadtxt(myfile, skiprows = 0)
        df = pd.DataFrame(data)
        df.columns = colArray
        return df
        
# Write out file as CSV for R load e.g. read.csv("my_R_file")
outfilename = sys.argv[2]
infilename = sys.argv[1]
print(infilename)
df = GSLIB2Dataframe(infilename)
df.to_csv(outfilename, sep=',',index=False)

