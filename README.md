# Separation of periodic orbits in the delay embedded space of chaotic attractors

This work explores the intersection of time-delay embeddings, periodic orbit theory, and
symbolic dynamics. Time-delay embeddings have been effectively applied to chaotic time-series
data, offering a principled method to reconstruct relevant information of the full attractor from
partial time-series observations. In this study, we investigate the structure of the unstable
periodic orbits of an attractor using time-delay embeddings. First, we embed time-series data
from a periodic orbit into a higher-dimensional space through the construction of Hankel matrix,
composed of time-shifted copies of the data. We then examine the influence of the width and
height of the Hankel matrix on the geometry of unstable periodic orbits in the delay-embedded
space. The right singular vectors of the Hankel matrix provide a basis for embedding the
periodic orbits. We observe that increasing the height of the Hankel matrix leads to a clear
separation of the periodic orbits into distinct clusters within the embedded space. Our analysis
characterizes these separated clusters and provides a mathematical framework to determine the
relative position of individual unstable periodic orbits in the embedded space. Additionally, we
present a modified formula to derive symbolic representation of distinct periodic orbits for a
specified sequence length, extending the Polyá-Redfield enumeration theorem.

The publication "Separation of periodic orbits in the delay embedded space of chaotic attractors" by Prerna Patil, Eurika Kaiser, JN Kutz and Steven Brunton is available on [arXiv](https://arxiv.org/abs/2411.13103)

# Code organization 
The data for unstable periodic orbits is available in the `Data/` folder. Two separate folders are available: 1.) `Lorenz/` 2.) `Rossler/`
The intial conditions for the UPOs can be obtained in these folders and can be used to reconstruct the UPOs. 

The folder `Codes/utils/` contains helper functions. These include construction of the Hankel matrix (`CreateHankelMatrix.m`), 
signal augmentation (`AugmentSignal.m`), create visualization of the attractor (`MakeAttractorPlot.m`) etc. 

# Lorenz attractor 
The separation of UPOs for the Lorenz attractor for the three cases: 1.) $`A^nB/AB^n`$ 2.) $`A^nB^n`$ 3.) UPOs with sequence length less than 8, 
can be obtained by the running the file `Codes/Lorenz/UnfoldingOfViswanathUPOs.m`. Selection of the cases can be done by commenting/uncommenting
the approprite lines between 35 to 37. The lenght of the time delay of the Hankel matrix can be set on line 129 to obtain different levels 
of separation. This code can be used to obtain Figure 6/8/10 (a-f) in the paper. 

Clustering of the UPOs based on the A-heavy vs. B-heavy orbits on the opposite ends of the unfolded attractor can be obtained by running 
the code `Codes/Lorenz/ColorCodedUPOs.m`. This code can be used to obtain Figure 9 in the paper. 

The comparison between UPOs based on their proximity in the embedded space vs. phase space can be obtained by running the code `Code/Lorenz/ProximityInSpaces.m`.
This code can be used to obtain Figure 7(a-b) in the publication.

# Rössler attractor
The separation of UPOs for the Rössler attractor for the two cases: 1.) $`A^nB/AB^n`$ 2.) UPOs with sequence length less than 8, can be
obtained by running the file `Codes/Rossler/ROSSLER_ParametersOfHankelMatrix.m`. This code can be used to obtain Figures 11/12 (a-f) 
in the paper. 