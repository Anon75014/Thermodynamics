# Thermodynamics of concentrated aqueous solutions
MATLAB code related to the submitted publication "Estimation of activity coefficients for aqueous organic redox-flow batteries: Theoretical basis and equations"

To reproduce figures that appear in the publication .pdf file, in order:

Figure 2, p.9:

In the folder "NaCl", please execute the files VisualisationSurfaceLPhi.m (excess enthalpy, subfigure a), VisualisationSurfaceJphi.m (excess heat capacity, subfigure b), VisualisationSurfacePhi.m (osmotic coefficient, subfigure c) and VisualisationSurfaceGamma.m (activity coefficient, subfigure d)

The input .txt files, ApparentHeatCapacityCpphiNaCl.txt, ExcessEnthalpyLPhiNaClMessikomer.txt, OsmoticCoefficientNaClGibbard.txt and ReferenceActivityNaCl_clean.txt, were generated manually from the fitting datasets, taken as reported from the peer-reviewed punlications cited in the manuscript.

Figure 4, p.11:

The figure is plotted using LaTeX tikz, calling the relevant .csv files in the NaCl folder, : phiNaClCalc.csv (lines) and phiNaClExp.csv (dots) for figure 4a, gammaNaClCalc.csv (lines) and gammaNaClExp.csv (dots) for figure 4b, and gammaNaClResiduals.csv (dots) and phiNaClResiduals.csv (dots) for figure C6.
The .csv files themselves are generated from GammaNaClValidation.m and PhiNaClValidation.m, respectively.

In turn, the matlab programs GammaNaClValidation.m and PhiNaClValidation.m use the text files ReferenceActivityNaCl_clean.txt and ReferenceOsmoticNaCl_clean_2.txt as input data. These text files were created from manual compilation of the fitting datasets, taken as reported from the peer-reviewed publications cited in the manuscript.

Figure 5, p.12:

Same as the previous figure, the .csv files were generated from VirialMatrix_KCl.m and VirialMatrix_CaCl2.m, in the KCl and CaCl2 folders, respectively

