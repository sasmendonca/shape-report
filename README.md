ROCS Shape Report Generator
This script generates a comprehensive PDF report based on ROCS (Rapid Overlay of Chemical Structures) calculations. It analyzes molecular shapes, colors, and 2D overlaps against a user-defined query, providing insights into structural similarities and distributions.

Features
Customizable Query: Define your own query molecule to tailor the analysis.

Color Field Integration: Utilize a color field file (.cff) to visualize the distribution of chemical properties.

2D Overlap Depiction: Generate 2D structure overlays and Tanimoto similarity scores for visual comparison.

Optimized for ROCS Data: Designed to process output from ROCS, leveraging user-defined score functions efficiently.

Prerequisites
OpenEye Python environment (oepythonenv)

ROCS (Rapid Overlay of Chemical Structures)

How to Execute
Activate your OpenEye Python environment:

Bash

conda activate oepython
Run the script with the following command:

Bash

python shapeoverlap2pdf.py -in <input_oeb.gz_file> -out <output_pdf_report> -cfffile <color_field_file.cff> -maxhits <number_of_hits> -depictsim

Parameters:

-in <input_oeb.gz_file>: Path to the gzipped OpenEye Binary (OEB) file containing ROCS results.

-out <output_pdf_report>: Desired name for the output PDF report.

-cfffile <color_field_file.cff>: Path to the color field file.

-maxhits <number_of_hits>: Maximum number of hits to include in the report.

-depictsim: Flag to enable depiction with 2D structure overlap and Tanimoto similarity.

Important Considerations for .cff File
The .cff file is crucial for accurate color representation. Pay close attention to the SMARTS codes within this file. If the SMARTS codes are too general, they may not accurately reflect the desired color distribution around your query molecule. It is highly recommended to modify the SMARTS codes to better represent the specific chemical features and their corresponding colors relevant to your analysis.
