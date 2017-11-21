# SedInConnect
Executable files are available under the "Release Tab" at the following link: https://github.com/HydrogeomorphologyTools/SedInConnect_2.3/releases

Open source program for the assessment of sediment connectivity as expressed in:

 - Crema, S., Cavalli, M., 2018. SedInConnect: a stand-alone, free and open source tool for the assessment of sediment connectivity. Comput. Geosci. 111, 39–45. https://doi.org/10.1016/j.cageo.2017.10.009
   link to the article: http://www.sciencedirect.com/science/article/pii/S0098300416308275

 - Cavalli, M., Trevisani, S., Comiti, F., Marchi, L., 2013. Geomorphometric assessment of spatial sediment connectivity in small Alpine catchments. Geomorphology, Sediment sources, source-to-sink fluxes and sedimentary budgets 188, 31–41. doi:10.1016/j.geomorph.2012.05.007
   link to the article: http://www.sciencedirect.com/science/article/pii/S0169555X12002267
   
 - Crema, S., Schenato, L., Goldin, B., Marchi, L., Cavalli, M., 2015. Toward the development of a stand-alone application for the assessment of sediment connectivity. Rendiconti Online Della Soc. Geol. Ital. 34, 58–61. doi:10.3301/ROL.2015.37
   link to the article: http://rendiconti.socgeol.it/244/fulltext.html?ida=1845


Contents:
Guidelines, source code and ancillary files.

Executables (created using Pyinstaller, https://github.com/pyinstaller/pyinstaller/wiki), are available under the "Release Tab" at the following link: https://github.com/HydrogeomorphologyTools/SedInConnect_2.3/releases


# SedInConnect Windows dependencies

To run SedInConnect you need to install TauDEM tools from: http://hydrology.usu.edu/taudem/taudem5/downloads2.html.

The complete Windows installer is fine. The installation will check and, in case, install GDAL and HPC Pack. Check your computer settings, especially Firewall and/or antivirus SW not to block the MPI/smpd calls.


# SedInConnect Linux dependencies

Linux release has been tested in Ubuntu 14 (LTS).

Dependencies are:
MPI (sudo apt-get install mpi), if the mpiexec command raises error one fix could be to install the dependencies (sudo apt-get install libgdk-pixbuf2.0-dev).

TauDEM functions.

To compile TauDEM tools refer to:
https://github.com/alexbruy/sextante-taudem/blob/master/README.TauDEM
working on an older version of TauDEM (5.0.6)