##Script to generate ROCS shape report based on color, shape, and 2d overlap by user defined query

to execute: (oepythonenv) python shapeoverlap2pdf.py -in oeb.gzfile -out pdfreport -cfffile colorfieldfile -maxhits xx -depictsim (depict-2d-overlap)

depictsim is a command to depict with 2d structure overlap and tanimoto similarity

Script optimized to generate shape report from ROCS with user defined query and score function used.

Pay attention at the .cff file, the SMARTS code is too general and should be modified to better represents the colors distribuition around the query.
