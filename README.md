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
The data for unstable periodic orbits is available in the 'Data/' folder. Two separate folders are available: 1.) 'Lorenz/' 2.) 'Rossler/'