# Thermodynamics of concentrated aqueous solutions
MATLAB code related to the submitted publication "Estimation of activity coefficients for aqueous organic redox-flow batteries: Theoretical basis and equations"

To reproduce figures that appear in the publication .pdf file, in order:

Figure 2, p.9:

In the folder "NaCl", please execute the files VisualisationSurfaceLPhi.m (excess enthalpy, subfigure a), VisualisationSurfaceJphi.m (excess heat capacity, subfigure b), VisualisationSurfacePhi.m (osmotic coefficient, subfigure c) and VisualisationSurfaceGamma.m (activity coefficient, subfigure d)

The relevant .txt files, which compile data from the literature (referenced in the paper) are called on line 78.

Figure 4, p.11:

The figure is plotted using LaTeX tikz, calling the relevant .csv files (gammaNaClCalc.csv, gammaNaClExp.csv, gammaNaClResiduals.csv, phiNaClCalc.csv, phiNaClExp.csv and phiNaClResiduals.csv).

In the NaCl folder, the .csv files are generated from GammaNaClValidation.m and PhiNaClValidation.m, respectively, which call the relevant .txt files which compile data from the literature (referenced in the paper)

Figure 5, p.12:

Same as the previous figure, the .csv files were generated from VirialMatrix_KCl.m and VirialMatrix_CaCl.m, in the KCl and CaCl2 folders

