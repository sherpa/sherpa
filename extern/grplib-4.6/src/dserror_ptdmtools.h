/*                                                                
**  Copyright (C) 1998,2007  Smithsonian Astrophysical Observatory 
*/                                                                

/*                                                                          */
/*  This program is free software; you can redistribute it and/or modify    */
/*  it under the terms of the GNU General Public License as published by    */
/*  the Free Software Foundation; either version 3 of the License, or       */
/*  (at your option) any later version.                                     */
/*                                                                          */
/*  This program is distributed in the hope that it will be useful,         */
/*  but WITHOUT ANY WARRANTY; without even the implied warranty of          */
/*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           */
/*  GNU General Public License for more details.                            */
/*                                                                          */
/*  You should have received a copy of the GNU General Public License along */
/*  with this program; if not, write to the Free Software Foundation, Inc., */
/*  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.             */
/*                                                                          */


#ifndef DSERROR_STRUCTS_H
#include "dserror_structs.h"
#endif
#define _DSERROR_PTDM_H

#define dsPTDMTOOLERROFFSET     -2000
#define dsPTDMTOOLNUMERRORS     96


/* dmhedit errors */
#define dsDMHEDITBADLINEERR (dsPTDMTOOLERROFFSET - 50 )
#define dsDMHEDITBADLINESEV dsERRSEVWARNING
#define dsDMHEDITBADLINESTDMSG "WARNING: The line '%s' in file '%s' is incorrectly formatted.\n"

#define dsDMHEDITNOKEYWRITEERR (dsPTDMTOOLERROFFSET - 51 )
#define dsDMHEDITNOKEYWRITESEV dsERRSEVWARNING
#define dsDMHEDITNOKEYWRITESTDMSG "WARNING: The key '%s' was not written.\n"

#define dsDMHEDITNOKEYADDERR (dsPTDMTOOLERROFFSET - 52 )
#define dsDMHEDITNOKEYADDSEV dsERRSEVWARNING
#define dsDMHEDITNOKEYADDSTDMSG "WARNING: The key '%s' was not added.\n"

#define dsDMHEDITKEYDELETEERR (dsPTDMTOOLERROFFSET - 53 )
#define dsDMHEDITKEYDELETESEV dsERRSEVWARNING
#define dsDMHEDITKEYDELETESTDMSG "WARNING: The key '%s' was not deleted.\n"

#define dsDMHEDITKEYMOVEAFTERERR (dsPTDMTOOLERROFFSET - 54 )
#define dsDMHEDITKEYMOVEAFTERSEV dsERRSEVWARNING
#define dsDMHEDITKEYMOVEAFTERSTDMSG "WARNING: Key '%s'  was not found, so moveafter failed.\n"

#define dsDMHEDITUNKNOWNOPPERR (dsPTDMTOOLERROFFSET - 55 )
#define dsDMHEDITUNKNOWNOPPSEV dsERRSEVWARNING
#define dsDMHEDITUNKNOWNOPPSTDMSG "WARNING: Operation '%s' unknown.\n"


/* dmsort errors */
#define dsDMSMISSTABLENAMEERR (dsPTDMTOOLERROFFSET - 200)
#define dsDMSMISSTABLENAMESEV dsERRSEVFATAL
#define dsDMSMISSTABLENAMESTDMSG "ERROR: Input sorting table name please!\n"

#define dsDMSINPUTFILEERR (dsPTDMTOOLERROFFSET - 201)
#define dsDMSINPUTFILESEV dsERRSEVFATAL
#define dsDMSINPUTFILESTDMSG "ERROR: Cannot sort image table %s.\n"


/* dmwritefef errors */
#define dsDMWFEFNOSCRPTFILNAMERR (dsPTDMTOOLERROFFSET - 250)
#define dsDMWFEFNOSCRPTFILNAMSEV  dsERRSEVFATAL
#define dsDMWFEFNOSCRPTFILNAMSTDMSG "ERROR: No file name given for script file.\n"

#define dsDMWFEFNOXFILNAMERR (dsPTDMTOOLERROFFSET - 251)
#define dsDMWFEFNOXFILNAMSEV  dsERRSEVFATAL
#define dsDMWFEFNOXFILNAMSTDMSG "ERROR: No file name given for axis file.\n"

#define dsDMWFEFNOLSTFILNAMERR (dsPTDMTOOLERROFFSET - 252)
#define dsDMWFEFNOLSTFILNAMSEV  dsERRSEVFATAL
#define dsDMWFEFNOLSTFILNAMSTDMSG "ERROR: No file name given for list file.\n"

#define dsDMWFEFFXNAMMISMATERR (dsPTDMTOOLERROFFSET - 253)
#define dsDMWFEFFXNAMMISMATSEV  dsERRSEVFATAL
#define dsDMWFEFFXNAMMISMATSTDMSG "ERROR: Axis name mismatch with %s.\n"

#define dsDMWFEFCANTSETXPROPERR (dsPTDMTOOLERROFFSET - 254)
#define dsDMWFEFCANTSETXPROPSEV  dsERRSEVFATAL
#define dsDMWFEFCANTSETXPROPSTDMSG "ERROR: Couldn't set axis name/min/max, units, or scale.\n"

#define dsDMWFEFXNOMISMATERR (dsPTDMTOOLERROFFSET - 255)
#define dsDMWFEFXNOMISMATSEV  dsERRSEVFATAL
#define dsDMWFEFXNOMISMATSTDMSG "ERROR: # of axes from ascfit (%d) larger than # of axes given (%d).\n"

#define dsDMWFEFCANTREADCONSTERR (dsPTDMTOOLERROFFSET - 256)
#define dsDMWFEFCANTREADCONSTSEV  dsERRSEVFATAL
#define dsDMWFEFCANTREADCONSTSTDMSG "ERROR: Couldn't read const parameters.\n"

#define dsDMWFEFMOREPARAMCOLNERR (dsPTDMTOOLERROFFSET - 257)
#define dsDMWFEFMOREPARAMCOLNSEV  dsERRSEVFATAL
#define dsDMWFEFMOREPARAMCOLNSTDMSG "ERROR: Too many non parameter columns in list file.\n"

#define dsDMWFEFXNOMISMATGVNERR (dsPTDMTOOLERROFFSET - 258)
#define dsDMWFEFXNOMISMATGVNSEV  dsERRSEVFATAL
#define dsDMWFEFXNOMISMATGVNSTDMSG "ERROR: # of axes in list file (%d) plus # of ascfit axes (%d) is different than number of axes given (%d).\n"

#define dsDMWFEFCANTSETPARAMERR (dsPTDMTOOLERROFFSET - 259)
#define dsDMWFEFCANTSETPARAMSEV  dsERRSEVFATAL
#define dsDMWFEFCANTSETPARAMSTDMSG "ERROR: Couldn't set parameter names/values.\n"

#define dsDMWFEFMISCNTPARAMERR (dsPTDMTOOLERROFFSET - 260)
#define dsDMWFEFMISCNTPARAMSEV  dsERRSEVFATAL
#define dsDMWFEFMISCNTPARAMSTDMSG "ERROR: Parameter miscount.\n"

#define dsDMWFEFCANTSETCOLNAMERR (dsPTDMTOOLERROFFSET - 261)
#define dsDMWFEFCANTSETCOLNAMSEV  dsERRSEVFATAL
#define dsDMWFEFCANTSETCOLNAMSTDMSG "ERROR: Couldn't set column names.\n"

#define dsDMWFEFFORKERRERR (dsPTDMTOOLERROFFSET - 262)
#define dsDMWFEFFORKERRSEV  dsERRSEVFATAL
#define dsDMWFEFFORKERRSTDMSG "ERROR: Fork error.\n"

/* ---------------------------------------- */
/* dmextract errors */
/* ---------------------------------------- */
/* ignored parameter-not in use */
#define dsDMEXTRACTIGNOREPARERR (dsPTDMTOOLERROFFSET - 350 )
#define dsDMEXTRACTIGNOREPARSEV dsERRSEVWARNING
#define dsDMEXTRACTIGNOREPARSTDMSG "WARNING: The parameter '%s' will be ignored.\n"
/* binning spec missing-not in use */
#define dsDMEXTRACTBINSPECMISSINGERR (dsPTDMTOOLERROFFSET - 351 )
#define dsDMEXTRACTBINSPECMISSINGSEV dsERRSEVFATAL
#define dsDMEXTRACTBINSPECMISSINGSTDMSG "ERROR: The binning specification is missing.\n"
/* wrong format of binning spec-not in use */
#define dsDMEXTRACTBINSPECMISSINGEQERR (dsPTDMTOOLERROFFSET - 352 )
#define dsDMEXTRACTBINSPECMISSINGEQSEV dsERRSEVFATAL
#define dsDMEXTRACTBINSPECMISSINGEQSTDMSG "ERROR: The binning specification is missing the '='.\n"
/*wrong format of binnning spec-not in use*/
#define dsDMEXTRACTBINPARMISSINGERR (dsPTDMTOOLERROFFSET - 353 )
#define dsDMEXTRACTBINPARMISSINGSEV dsERRSEVWARNING
#define dsDMEXTRACTBINPARMISSINGSTDMSG "WARNING: The binning parameter '%s' is not specified by the user.\n"
/*falure to find the bin-not in use*/
#define dsDMEXTRACTFINDBINPARERR (dsPTDMTOOLERROFFSET - 354 )
#define dsDMEXTRACTFINDBINPARSEV dsERRSEVFATAL
#define dsDMEXTRACTFINDBINPARSTDMSG "ERROR: A value the binning parameter '%s' was not found.\n"
/*bins are backwords*/
#define dsDMEXTRACTMINGTMAXERR (dsPTDMTOOLERROFFSET - 355 )
#define dsDMEXTRACTMINGTMAXSEV dsERRSEVFATAL
#define dsDMEXTRACTMINGTMAXSTDMSG "ERROR: Minbin is greater than maxbin. Check parameters.\n"
/*binning is 0-not in use*/
#define dsDMEXTRACTBINSIZENEGERR (dsPTDMTOOLERROFFSET - 356 )
#define dsDMEXTRACTBINSIZENEGSEV dsERRSEVFATAL
#define dsDMEXTRACTBINSIZENEGSTDMSG "ERROR: Bin factor or array size is <= 0.\n"
/*bad option*/
#define dsDMEXTRACTBADOPTERR (dsPTDMTOOLERROFFSET - 357)
#define dsDMEXTRACTBADOPTSEV dsERRSEVFATAL
#define dsDMEXTRACTBADOPTSTDMSG "ERROR: Unrecognized option: %s\n"

/*background error */
#define dsDMEXTRACTBKGNUMERR (dsPTDMTOOLERROFFSET - 358)
#define dsDMEXTRACTBKGNUMSEV dsERRSEVFATAL
#define dsDMEXTRACTBKGNUMSTDMSG "ERROR: Background regions must have the same number as the input region\n\tunless only 1 background region is supplied.\n"

/*data file error */
#define dsDMEXTRACTDATAREADERR (dsPTDMTOOLERROFFSET - 359)
#define dsDMEXTRACTDATAREADSEV dsERRSEVFATAL
#define dsDMEXTRACTDATAREADSTDMSG "ERROR: Problem reading the data from %s.\n"

/*process all files */
#define dsDMEXTRACTPROFILEERR (dsPTDMTOOLERROFFSET - 360)
#define dsDMEXTRACTPROFILESEV dsERRSEVFATAL
#define dsDMEXTRACTPROFILESTDMSG "ERROR: Failed to process some files.\n"

/*region parse */
#define dsDMEXTRACTREGPARSEERR (dsPTDMTOOLERROFFSET - 361)
#define dsDMEXTRACTREGPARSESEV dsERRSEVFATAL
#define dsDMEXTRACTREGPARSESTDMSG "ERROR: Failed to parse the supplied regions. Please check the format.\n"

/*step size negative */
#define dsDMEXTRACTSTEPNEGERR (dsPTDMTOOLERROFFSET -362)
#define dsDMEXTRACTSTEPNEGSEV dsERRSEVWARNING
#define dsDMEXTRACTSTEPNEGSTDMSG "WARNING: Step size negative: sign ignored\n"

/*step size zero */
#define dsDMEXTRACTSTEPZEROERR (dsPTDMTOOLERROFFSET - 363)
#define dsDMEXTRACTSTEPZEROSEV dsERRSEVWARNING
#define dsDMEXTRACTSTEPZEROSTDMSG "WARNING: Step size zero; set to +1.0\n"

/*area warning */
#define dsDMEXTRACTNOSKYWERR (dsPTDMTOOLERROFFSET - 364)
#define dsDMEXTRACTNOSKYWSEV dsERRSEVWARNING
#define dsDMEXTRACTNOSKYWSTDMSG "WARNING: No SKY or POS(X,Y) column: could not calculate area\n"

/* wcs reg warning */
#define dsDMEXTRACTREGWCSWERR (dsPTDMTOOLERROFFSET - 365)
#define dsDMEXTRACTREGWCSWSEV dsERRSEVWARNING
#define dsDMEXTRACTREGWCSWSTDMSG "WARNING: Expected 2 dimensions for the sky descriptor. Skipping translation of region to WCS.\n"

/*header key problem */
#define dsDMEXTRACTLOADKEYWERR (dsPTDMTOOLERROFFSET - 366)
#define dsDMEXTRACTLOADKEYWSEV dsERRSEVWARNING
#define dsDMEXTRACTLOADKEYWSTDMSG "WARNING: Problems loading header keyword %s into output file. Skipping keyword...\n"

/*default Extraction */
#define dsDMEXTRACTDEFAULTEXWERR (dsPTDMTOOLERROFFSET - 367)
#define dsDMEXTRACTDEFAULTEXWSEV dsERRSEVWARNING
#define dsDMEXTRACTDEFAULTEXWSTDMSG "WARNING: No extraction given for file %s, assuming default PI extraction\n"

/*header read error */
#define dsDMEXTRACTHEADERREADERR (dsPTDMTOOLERROFFSET - 368)
#define dsDMEXTRACTHEADERREADSEV dsERRSEVFATAL
#define dsDMEXTRACTHEADERREADSTDMSG "ERROR: Error reading header in file %s.\n"

/*exposure stack error */
#define dsDMEXTRACTEXPSTACKERR (dsPTDMTOOLERROFFSET - 369)
#define dsDMEXTRACTEXPSTACKSEV dsERRSEVFATAL
#define dsDMEXTRACTEXPSTACKSTDMSG "ERROR: The number of exposure map files and input files do not match.\n"

/*background stack error */	
#define dsDMEXTRACTBKGSTACKERR (dsPTDMTOOLERROFFSET - 370)
#define dsDMEXTRACTBKGSTACKSEV dsERRSEVFATAL
#define dsDMEXTRACTBKGSTACKSTDMSG "ERROR: The number of background files and input files do not match.\n"

/*background exp stack error */	
#define dsDMEXTRACTBKGEXPSTACKERR (dsPTDMTOOLERROFFSET - 371)
#define dsDMEXTRACTBKGEXPSTACKSEV dsERRSEVFATAL
#define dsDMEXTRACTBKGEXPSTACKSTDMSG "ERROR: The number of background files and background exposure map files do not match.\n"

/*more than 1 component warning */
#define dsDMEXTRACTREGCOMPWERR (dsPTDMTOOLERROFFSET - 372)
#define dsDMEXTRACTREGCOMPWSEV dsERRSEVWARNING
#define dsDMEXTRACTREGCOMPWSTDMSG "WARNING:Region #%d contains more than 1 component. Only the first component will be described in the region columns of the output file.\n"

/* no regions on source file */
#define dsDMEXTRACTNOREGERR (dsPTDMTOOLERROFFSET - 373)
#define dsDMEXTRACTNOREGSEV dsERRSEVFATAL
#define dsDMEXTRACTNOREGSTDMSG "ERROR: Source image has no extraction regions supplied.\n"

/*variance error */
#define dsDMEXTRACTVARIMGERR (dsPTDMTOOLERROFFSET - 374)
#define dsDMEXTRACTVARIMGSEV dsERRSEVWARNING
#define dsDMEXTRACTVARIMGSTDMSG "WARNING: Can't open variance file \"%s\". Using gaussian errors.\n"

/*background no regions */
#define dsDMEXTRACTBKGNOREGERR (dsPTDMTOOLERROFFSET - 375)
#define dsDMEXTRACTBKGNOREGSEV dsERRSEVWARNING
#define dsDMEXTRACTBKGNOREGSTDMSG "WARNING: Background data NOT extracted over a region. Setting background to 0.0.\n"

/*subarray error */
#define dsDMEXTRACTSUBARRAYERR (dsPTDMTOOLERROFFSET - 376)
#define dsDMEXTRACTSUBARRAYSEV dsERRSEVFATAL
#define dsDMEXTRACTSUBARRAYSTDMSG "ERROR: Failed to read in subarray\n"

/* norm zero warning */
#define dsDMEXTRACTNOBKGERR (dsPTDMTOOLERROFFSET - 377)
#define dsDMEXTRACTNOBKGSEV dsERRSEVWARNING
#define dsDMEXTRACTNOBKGSTDMSG "WARNING:Setting background normalization to 0 will force the background to be 0.\n"

/*variance error */
#define dsDMEXTRACTVARPHAERR (dsPTDMTOOLERROFFSET - 378)
#define dsDMEXTRACTVARPHASEV dsERRSEVWARNING
#define dsDMEXTRACTVARPHASTDMSG "WARNING: Variance files are not yet supported for histogram/pha extractions.\n\tUsing gaussian errors.\n"
/* ---------------------------------------- */
/* dmregrid errors */
#define dsDMREGRIDINCOMPLETEBINSPECERR (dsPTDMTOOLERROFFSET - 450 )
#define dsDMREGRIDINCOMPLETEBINSPECSEV dsERRSEVFATAL
#define dsDMREGRIDINCOMPLETEBINSPECSTDMSG  "ERROR: Binning specification is incomplete (missing a colon, comma, or quote.\n"

#define dsDMREGRIDBADNUMBINSPECERR (dsPTDMTOOLERROFFSET - 451 )
#define dsDMREGRIDBADNUMBINSPECSEV dsERRSEVFATAL
#define dsDMREGRIDBADNUMBINSPECSTDMSG  "ERROR: The number of binning specifications is incompatible with the number of input files.\n"

#define dsDMREGRIDIMAGEEXTENTDIFFERERR (dsPTDMTOOLERROFFSET - 452 )
#define dsDMREGRIDIMAGEEXTENTDIFFERSEV dsERRSEVWARNING
#define dsDMREGRIDIMAGEEXTENTDIFFERSTDMSG "WARNING: The extent of image %d differs from the extent of  thefirst image.\n"

#define dsDMREGRIDMINGTMAXERR (dsPTDMTOOLERROFFSET - 453 )
#define dsDMREGRIDMINGTMAXSEV dsERRSEVFATAL
#define dsDMREGRIDMINGTMAXSTDMSG "ERROR: Minbin is greater than maxbin. Check parameters.\n"

#define dsDMREGRIDBINSIZENEGERR (dsPTDMTOOLERROFFSET - 454 )
#define dsDMREGRIDBINSIZENEGSEV dsERRSEVFATAL
#define dsDMREGRIDBINSIZENEGSTDMSG "ERROR: Bin factor or array size is <= 0.\n"


/* dmgroup errors */
#define dsDMGROUPUNSUPPGTYPEERR  ( dsPTDMTOOLERROFFSET - 500 )
#define dsDMGROUPUNSUPPGTYPESEV dsERRSEVFATAL
#define dsDMGROUPUNSUPPGTYPESTDMSG "ERROR: input group type is not supported.\n"

#define dsDMGROUPMINBINGTMAXBINERR ( dsPTDMTOOLERROFFSET - 501 )
#define dsDMGROUPMINBINGTMAXBINSEV dsERRSEVFATAL
#define dsDMGROUPMINBINGTMAXBINSTDMSG "ERROR: Min bin is greater than max bin. Check parameters.\n"

#define dsDMGROUPSTEPNEGERR ( dsPTDMTOOLERROFFSET - 502 )  
#define dsDMGROUPSTEPNEGSEV dsERRSEVWARNING
#define dsDMGROUPSTEPNEGSTDMSG "WARNING: The bin size is negative, so changing sign.\n"

#define dsDMGROUPSTEPZEROERR (dsPTDMTOOLERROFFSET - 503 )  
#define dsDMGROUPSTEPZEROSEV  dsERRSEVWARNING
#define dsDMGROUPSTEPZEROSTDMSG "WARNING: The given bin size is zero; it has been reset to %f.\n"

#define dsDMGROUPOVERLAPBINSPECERR ( dsPTDMTOOLERROFFSET - 504 )
#define dsDMGROUPOVERLAPBINSPECSEV dsERRSEVFATAL
#define dsDMGROUPOVERLAPBINSPECSTDMSG "ERROR: At least 2 binning specifications overlap. Check binning specifications.\n"

#define dsDMGROUPMINBINERR ( dsPTDMTOOLERROFFSET - 505 )
#define dsDMGROUPMINBINSEV dsERRSEVWARNING
#define dsDMGROUPMINBINSTDMSG "ERROR: The minimum binning channel is out of range; it has been reset to %d.\n"

#define dsDMGROUPMAXBINERR ( dsPTDMTOOLERROFFSET - 506 )
#define dsDMGROUPMAXBINSEV dsERRSEVWARNING
#define dsDMGROUPMAXBINSTDMSG "ERROR: The maximum binning channel is out of range; it has been reset to %d.\n"

#define dsDMGROUPNOBINSPECERR ( dsPTDMTOOLERROFFSET - 507 )
#define dsDMGROUPNOBINSPECSEV dsERRSEVFATAL
#define dsDMGROUPNOBINSPECSTDMSG "ERROR: No binning specification was given.\n"

#define dsDMGROUPBINSPECMAXERR ( dsPTDMTOOLERROFFSET - 508 )
#define dsDMGROUPBINSPECMAXSEV dsERRSEVWARNING
#define dsDMGROUPBINSPECMAXSTDMSG "WARNING: The maximum channel given in the binning specification ( %f ) is too large for the \n given bin size ( %f ) and value in last channel of input file ( %d ); resetting the maximum \n channel to %f.\n" 

#define dsDMGROUPSTEPMAXERR ( dsPTDMTOOLERROFFSET - 509 )
#define dsDMGROUPSTEPMAXSEV dsERRSEVWARNING
#define dsDMGROUPSTEPMAXSTDMSG "WARNING: The bin size given in the binning specification is too great; resetting the bin size to %f.\n"

#define dsDMGROUPBADPARAMERR ( dsPTDMTOOLERROFFSET - 510 )
#define dsDMGROUPBADPARAMSEV dsERRSEVFATAL
#define dsDMGROUPBADPARAMSTDMSG "ERROR: At least one input parameter has an invalid value.\n"

#define dsDMGROUPBADDATAORDERERR ( dsPTDMTOOLERROFFSET - 511 )
#define dsDMGROUPBADDATAORDERSEV dsERRSEVFATAL
#define dsDMGROUPBADDATAORDERSTDMSG "ERROR: Data column is not monotonically increasing.\n" 

#define dsDMGROUPZEROWIDTHERR ( dsPTDMTOOLERROFFSET - 512 )
#define dsDMGROUPZEROWIDTHSEV dsERRSEVWARNING
#define dsDMGROUPZEROWIDTHSTDMSG "WARNING: The calculated bin width rounds to zero.\nIt will be reset to 1.\n"

#define dsDMGROUPINVALIDBINERR ( dsPTDMTOOLERROFFSET - 513 )
#define dsDMGROUPINVALIDBINSEV dsERRSEVFATAL
#define dsDMGROUPINVALIDBINSTDMSG "ERROR: A bin boundary is invalid.\nMake sure the binspec fits within the bounds of the data.\n"

#define dsDMGROUPLOWERBOUNDERR ( dsPTDMTOOLERROFFSET - 514 )
#define dsDMGROUPLOWERBOUNDSEV dsERRSEVFATAL
#define dsDMGROUPLOWERBOUNDSTDMSG "ERROR: grp_priv.c:lower_bound(): No data greater than or equal to given value.\n"

#define dsDMGROUPUPPERBOUNDERR ( dsPTDMTOOLERROFFSET - 515 )
#define dsDMGROUPUPPERBOUNDSEV dsERRSEVFATAL
#define dsDMGROUPUPPERBOUNDSTDMSG "ERROR: grp_priv.c:upper_bound(): No data less than or equal to given value.\n"

#define dsDMGROUPBINONROWSERR ( dsPTDMTOOLERROFFSET - 516 )
#define dsDMGROUPBINONROWSSEV dsERRSEVWARNING
#define dsDMGROUPBINONROWSSTDMSG "WARNING: The %s parameter is empty; will proceed to bin on rows.\n"

#define dsDMGROUPMISSINGPARAMERR ( dsPTDMTOOLERROFFSET - 517 )
#define dsDMGROUPMISSINGPARAMSEV dsERRSEVFATAL
#define dsDMGROUPMISSINGPARAMSTDMSG "ERROR: Grouping method requires parameter %s, which is missing or out of bounds.\n"

#define dsDMGROUPEXTRAGROUPSERR ( dsPTDMTOOLERROFFSET - 518 )
#define dsDMGROUPEXTRAGROUPSSEV dsERRSEVWARNING
#define dsDMGROUPEXTRAGROUPSSTDMSG "WARNING: More groups produced than requested.\n"

#define dsDMGROUPTOOFEWGROUPSERR ( dsPTDMTOOLERROFFSET - 519 )
#define dsDMGROUPTOOFEWGROUPSSEV dsERRSEVFATAL
#define dsDMGROUPTOOFEWGROUPSSTDMSG "ERROR: Fewer groups produced than requested.\n"

#define dsDMGROUPZEROERRORERR ( dsPTDMTOOLERROFFSET - 520 )
#define dsDMGROUPZEROERRORSEV dsERRSEVWARNING
#define dsDMGROUPZEROERRORSTDMSG "WARNING: The supplied error column contains zero-valued data."


/* dmtype2split errors */
#define dsDMTYPE2SPLITBADNUMOFROWSERR ( dsPTDMTOOLERROFFSET - 550 )
#define dsDMTYPE2SPLITBADNUMOFROWSSEV    dsERRSEVFATAL
#define dsDMTYPE2SPLITBADNUMOFROWSSTDMSG  "ERROR:  The number of rows does not equal the number of output files .\n"

#define dsDMTYPE2SPLITBADROWNUMERR ( dsPTDMTOOLERROFFSET - 551 )
#define dsDMTYPE2SPLITBADROWNUMSEV    dsERRSEVWARNING
#define dsDMTYPE2SPLITBADROWNUMSTDMSG  "ERROR:  A requested row has an invalid numerical value, proceeding to next row.\n"

#define dsDMTYPE2SPLITROWNANERR ( dsPTDMTOOLERROFFSET - 552 )
#define dsDMTYPE2SPLITROWNANSEV    dsERRSEVFATAL
#define dsDMTYPE2SPLITROWNANSTDMSG  "ERROR:  A requested row has an non-numerical value.\n"

#define dsDMTYPE2SPLITBADROWERR ( dsPTDMTOOLERROFFSET - 553 )
#define dsDMTYPE2SPLITBADROWSEV    dsERRSEVFATAL
#define dsDMTYPE2SPLITBADROWSTDMSG  "ERROR:  A requested row has an value which could not be interpreted.\n"


/* DMDIFF ERRORS */
#define dsDMDIFFLINEDIFFERR ( dsPTDMTOOLERROFFSET - 600 )
#define dsDMDIFFLINEDIFFSEV    dsERRSEVWARNING
#define dsDMDIFFLINEDIFFSTDMSG  "\n WARNING: %ld line(s) different.\n\n"

#define dsDMDIFFUNITDIFFERR ( dsPTDMTOOLERROFFSET - 601 )
#define dsDMDIFFUNITDIFFSEV    dsERRSEVWARNING
#define dsDMDIFFUNITDIFFSTDMSG  "WARNING: Units are NOT the same in both files!\n"

#define dsDMDIFFUNITNOTEQERR ( dsPTDMTOOLERROFFSET - 602 )
#define dsDMDIFFUNITNOTEQSEV    dsERRSEVWARNING
#define dsDMDIFFUNITNOTEQSTDMSG  "WARNING: File %d UNIT does not equal \"%s\"\n         Refer to tolfile for expected value.\n"

#define dsDMDIFFCOMMENTDIFFERR ( dsPTDMTOOLERROFFSET - 603 )
#define dsDMDIFFCOMMENTDIFFSEV    dsERRSEVWARNING
#define dsDMDIFFCOMMENTDIFFSTDMSG  "WARNING: Comments are NOT the same in both files!\n"

#define dsDMDIFFCOMMENTNOTEQERR ( dsPTDMTOOLERROFFSET - 604 )
#define dsDMDIFFCOMMENTNOTEQSEV    dsERRSEVWARNING
#define dsDMDIFFCOMMENTNOTEQSTDMSG  "WARNING: File %d COMMENT does not equal \"%s\"\n         Refer to tolfile for expected value.\n"

#define dsDMDIFFLOWPARAMSERR ( dsPTDMTOOLERROFFSET - 605 )
#define dsDMDIFFLOWPARAMSSEV    dsERRSEVWARNING
#define dsDMDIFFLOWPARAMSSTDMSG  "WARNING: Not enough parameters in line %d of the tolfile.\n\n"

#define dsDMDIFFDATATYPEDIFFERR ( dsPTDMTOOLERROFFSET - 606 )
#define dsDMDIFFDATATYPEDIFFSEV    dsERRSEVWARNING
#define dsDMDIFFDATATYPEDIFFSTDMSG  "WARNING: keyword %s datatypes do not match in each file.\n"

#define dsDMDIFFVECTORDIFFERR ( dsPTDMTOOLERROFFSET - 607 )
#define dsDMDIFFVECTORDIFFSEV    dsERRSEVWARNING
#define dsDMDIFFVECTORDIFFSTDMSG  "WARNING: Vectors are DIFFERENT SIZES!\n"

#define dsDMDIFFROWSDIFFERR ( dsPTDMTOOLERROFFSET - 608 )
#define dsDMDIFFROWSDIFFSEV    dsERRSEVWARNING
#define dsDMDIFFROWSDIFFSTDMSG  " WARNING: Number of rows are different ... skipping."

#define dsDMDIFFSSVALSDIFFERR ( dsPTDMTOOLERROFFSET - 609 )
#define dsDMDIFFSSVALSDIFFSEV    dsERRSEVWARNING
#define dsDMDIFFSSVALSDIFFSTDMSG  "       WARNING:  num of subspace vals are different in each file\n"

#define dsDMDIFFVALSDIFFERR ( dsPTDMTOOLERROFFSET - 610 )
#define dsDMDIFFVALSDIFFSEV    dsERRSEVWARNING
#define dsDMDIFFVALSDIFFSTDMSG  "       WARNING: values are different in row:  %ld\n"

#define dsDMDIFFARRAYSIZEDIFFERR ( dsPTDMTOOLERROFFSET - 611 )
#define dsDMDIFFARRAYSIZEDIFFSEV    dsERRSEVWARNING
#define dsDMDIFFARRAYSIZEDIFFSTDMSG  "WARNING: Array sizes are different ... skipping.\n"

#define dsDMDIFFDATADIFFERR ( dsPTDMTOOLERROFFSET - 612 )
#define dsDMDIFFDATADIFFSEV    dsERRSEVWARNING
#define dsDMDIFFDATADIFFSTDMSG  "WARNING: data do not match in line %d. %d != %d\n"

#define dsDMDIFFDATANOTEQERR ( dsPTDMTOOLERROFFSET - 613 )
#define dsDMDIFFDATANOTEQSEV    dsERRSEVWARNING
#define dsDMDIFFDATANOTEQSTDMSG  "WARNING: data do not match. %d != %d\n"

#define dsDMDIFFSSDATATYPEDIFFERR ( dsPTDMTOOLERROFFSET - 614 )
#define dsDMDIFFSSDATATYPEDIFFSEV    dsERRSEVWARNING
#define dsDMDIFFSSDATATYPEDIFFSTDMSG  "WARNING: subspace column datatypes do not match.\n"

#define dsDMDIFFBADNVALSERR ( dsPTDMTOOLERROFFSET - 615 )
#define dsDMDIFFBADNVALSSEV    dsERRSEVFATAL
#define dsDMDIFFBADNVALSSTDMSG  "       ERROR: problem w/ nvalues: file1 (%ld) file2 (%ld)\n"
