########################################
# README file for SDHCALAnalysis package
# @author Eté Rémi
# @version 1.0.0 18/02/2015
# @copyright IPNL, CNRS
########################################


//SDHCALAnalysis


This file is part of SDHCALAnalysis libraries.

SDHCALAnalysis is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

SDHCALAnalysis is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with SDHCALAnalysis.  If not, see <http://www.gnu.org/licenses/>.


OVERVIEW:
==========

  The package consists mainly in :
  
    * ROOT macros that you can find in the 'macros' directory
    * Marlin processors that you can load with Marlin by exporting the MARLIN_DLL env variable the package library libSDHCALAnalysis.
    
    Simply do :
      $ export MARLIN_DLL=path/to/libSDHCALAnalysis.so
      $ Marlin path/to/MarlinXmlFile.xml
    
    The corresponding marlin xml file can be found in the 'xml' directory

INSTALL:
=========

See INSTALL file in the main directory
for more informations about installation.


Please, send comments or report bug via the github interface :

    https://github.com/SDHCAL/SDHCALAnalysis/issues
