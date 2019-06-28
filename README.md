# SedInConnect
Executable files are available under the "Release Tab" at the following link: https://github.com/HydrogeomorphologyTools/SedInConnect_2.3/releases

Open source program for the assessment of sediment connectivity as expressed in:

 - Crema, S., Cavalli, M., 2018. SedInConnect: a stand-alone, free and open source tool for the assessment of sediment connectivity. Comput. Geosci. 111, 39–45. https://doi.org/10.1016/j.cageo.2017.10.009

 - Cavalli, M., Trevisani, S., Comiti, F., Marchi, L., 2013. Geomorphometric assessment of spatial sediment connectivity in small Alpine catchments. Geomorphology, Sediment sources, source-to-sink fluxes and sedimentary budgets 188, 31–41. doi:10.1016/j.geomorph.2012.05.007
	http://www.sciencedirect.com/science/article/pii/S0169555X12002267
   
 - Crema, S., Schenato, L., Goldin, B., Marchi, L., Cavalli, M., 2015. Toward the development of a stand-alone application for the assessment of sediment connectivity. Rendiconti Online Della Soc. Geol. Ital. 34, 58–61. doi:10.3301/ROL.2015.37
	http://rendiconti.socgeol.it/244/fulltext.html?ida=1845


Contents:
Guidelines, source code and ancillary files.

Executables (created using Pyinstaller, https://github.com/pyinstaller/pyinstaller/wiki), are available under the "Release Tab" at the following link: https://github.com/HydrogeomorphologyTools/SedInConnect_2.3/releases


# SedInConnect Windows dependencies

To run SedInConnect you need to install TauDEM tools from: http://hydrology.usu.edu/taudem/taudem5/downloads2.html.

The complete Windows installer is fine. The installation will check and, in case, install GDAL and HPC Pack. Check your computer settings, especially Firewall and/or antivirus SW not to block the MPI/smpd calls.


# SedInConnect Linux dependencies

SedInConnect_2.3 Ubuntu16 (LTS) 64-bit testing release
Tested on Ubuntu 16.04.6 xenial

Dependencies: install TauDEM and mpi and ensure they work invoked by the command line (avoid calling taudem and mpiexec with sudo privilegies). Follow this link https://unix.stackexchange.com/questions/346775/how-to-solve-building-error-where-a-variable-was-not-declared-in-this-scope and/or run the following shell script https://sites.google.com/site/geoluislopes/taudem_ubuntu.tar.bz2?attredirects=0&d=1 to install both TauDEM and mpi

Uncompress the herein provided archive into a unique folder, change execute permissions to all the files in the folder by typing "sudo chmod -R a+rwx /path/to/folder", then run "./SedInConnect_2_3_lx_64" without sudo privilegies.

Have fun!
