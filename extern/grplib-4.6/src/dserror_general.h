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
#define _DSERROR_GENERAL_H

#define dsGENERALERROFFSET     0
#define dsGENERALNUMERRORS     92

#define dsGENCOLKEYOFFSET (dsGENERALERROFFSET - 25)
#define dsGENFIOOFFSET (dsGENCOLKEYOFFSET - 25)
#define dsGENDMOFFSET (dsGENFIOOFFSET - 25)
#define dsGENPARAMOFFSET (dsGENDMOFFSET - 50)
#define dsGENASCFITOFFSET (dsGENPARAMOFFSET - 25)
#define dsGENTIMEOFFSET (dsGENASCFITOFFSET - 25)
#define dsGENSTACKOFFSET (dsGENTIMEOFFSET - 25)
#define dsGENCALDBOFFSET (dsGENSTACKOFFSET - 25)

/* Miscellaneous errors */

#define dsNOERR dsGENERALERROFFSET
#define dsNOERRSEV  dsERRSEVNONE
#define dsNOERRSTDMSG ""

#define dsGENERICERR (dsGENERALERROFFSET - 1)
#define dsGENERICSEV dsERRSEVWARNING
#define dsGENERICSTDMSG "WARNING: An unspecified error has occurred.\n"

#define dsUNDEFERR (dsGENERALERROFFSET - 2)
#define dsUNDEFSEV dsERRSEVWARNING
#define dsUNDEFSTDMSG "WARNING: The error code '%s' is not defined.\n"

#define dsERRNOTFOUNDERR (dsGENERALERROFFSET - 3)
#define dsERRNOTFOUNDSEV dsERRSEVWARNING
#define dsERRNOTFOUNDSTDMSG "WARNING: Specified error was not found in the error list.\n"

#define dsDIVIDEBYZEROERR (dsGENERALERROFFSET - 4)
#define dsDIVIDEBYZEROSEV dsERRSEVFATAL
#define dsDIVIDEBYZEROSTDMSG "ERROR: Cannot divide by zero! \n"

#define dsALLOCERR (dsGENERALERROFFSET - 5)
#define dsALLOCSEV dsERRSEVFATAL
#define dsALLOCSTDMSG "ERROR: Could not allocate memory.\n"

#define dsALLOC2ERR (dsGENERALERROFFSET - 6)
#define dsALLOC2SEV dsERRSEVFATAL
#define dsALLOC2STDMSG "ERROR: Could not allocate memory for '%s'.\n"

#define dsALLOCOBJERR (dsGENERALERROFFSET - 7)
#define dsALLOCOBJSEV dsERRSEVFATAL
#define dsALLOCOBJSTDMSG "ERROR: Could not allocate new '%s' object.\n"

#define dsINITLIBERR (dsGENERALERROFFSET - 8)
#define dsINITLIBSEV dsERRSEVFATAL
#define dsINITLIBSTDMSG "ERROR: Could not initialize '%s' library.\n"

#define dsNULLPTRERR (dsGENERALERROFFSET - 9)
#define dsNULLPTRSEV dsERRSEVWARNING
#define dsNULLPTRSTDMSG "WARNING: A null ptr was received.\n"

#define dsNULLSTRINGERR (dsGENERALERROFFSET - 10)
#define dsNULLSTRINGSEV dsERRSEVWARNING
#define dsNULLSTRINGSTDMSG "WARNING: A null string was received.\n"

#define dsNOTIMPLMERR (dsGENERALERROFFSET - 11)
#define dsNOTIMPLMSEV dsERRSEVWARNING
#define dsNOTIMPLMSTDMSG "WARNING: The following functionality is not currently implemented:\n\t'%s'\n"

#define dsBOUNDSERR (dsGENERALERROFFSET - 12)
#define dsBOUNDSSEV dsERRSEVFATAL
#define dsBOUNDSSTDMSG "ERROR: Array '%s' bounds exceeded"

#define dsAUTONAMEERR (dsGENERALERROFFSET - 13)
#define dsAUTONAMESEV dsERRSEVFATAL
#define dsAUTONAMESTDMSG "ERROR: This tool cannot use autonaming. Please use an output file name.\n"

#define dsOBJSELFASSIGNERR (dsGENERALERROFFSET - 14)
#define dsOBJSELFASSIGNSEV dsERRSEVFATAL
#define dsOBJSELFASSIGNSTDMSG "ERROR: Attempting self assignment of a '%s' object\n"


/* Column and keyword errors in data files - mostly FITS, StSCI tables, or QPOES */
#define dsUNSUPORTTYPEERR (dsGENCOLKEYOFFSET - 1)
#define dsUNSUPORTTYPESEV dsERRSEVFATAL
#define dsUNSUPORTTYPESTDMSG "ERROR: An unsupported datatype was found in file '%s'.\n"

#define dsCOLUMNTYPEERR (dsGENCOLKEYOFFSET - 2)
#define dsCOLUMNTYPESEV dsERRSEVFATAL
#define dsCOLUMNTYPESTDMSG "ERROR: The column '%s' in file '%s' should have the datatype '%s'.\n"

#define dsKEYWORDTYPEERR (dsGENCOLKEYOFFSET - 3)
#define dsKEYWORDTYPESEV dsERRSEVFATAL
#define dsKEYWORDTYPESTDMSG "ERROR: The keyword '%s' in file '%s' should have the datatype '%s'.\n"

#define dsFINDCOLUMNERR (dsGENCOLKEYOFFSET - 4)
#define dsFINDCOLUMNSEV dsERRSEVFATAL
#define dsFINDCOLUMNSTDMSG "ERROR: Column '%s' was not found in file '%s'.\n"

#define dsFINDKEYWORDFERR (dsGENCOLKEYOFFSET - 5)
#define dsFINDKEYWORDFSEV dsERRSEVFATAL
#define dsFINDKEYWORDFSTDMSG "ERROR: Keyword '%s' was not found in file '%s'.\n"

#define dsFINDKEYWORDWERR (dsGENCOLKEYOFFSET - 6)
#define dsFINDKEYWORDWSEV dsERRSEVWARNING
#define dsFINDKEYWORDWSTDMSG "WARNING: Keyword '%s' was not found in file '%s', continuing without it.\n"

#define dsCREATECOLUMNERR (dsGENCOLKEYOFFSET - 7)
#define dsCREATECOLUMNSEV dsERRSEVFATAL
#define dsCREATECOLUMNSTDMSG "ERROR: Could not create column '%s' in file '%s'.\n"

#define dsCREATEKEYWORDERR (dsGENCOLKEYOFFSET - 8)
#define dsCREATEKEYWORDSEV dsERRSEVFATAL
#define dsCREATEKEYWORDSTDMSG "ERROR: Could not create keyword '%s' in file '%s'.\n"


/* file I/O errors */

#define dsOPENFILEFERR (dsGENFIOOFFSET - 1)
#define dsOPENFILEFSEV dsERRSEVFATAL
#define dsOPENFILEFSTDMSG "ERROR: Could not open file '%s'.\n"

#define dsOPENFILEWERR (dsGENFIOOFFSET - 2)
#define dsOPENFILEWSEV dsERRSEVWARNING
#define dsOPENFILEWSTDMSG "WARNING: Could not open file '%s', continuing without it.\n"

#define dsINPUTEXISTSFERR (dsGENFIOOFFSET - 3)
#define dsINPUTEXISTSFSEV dsERRSEVFATAL
#define dsINPUTEXISTSFSTDMSG "ERROR: Input file '%s' does not exist.\n"

#define dsINPUTEXISTSWERR (dsGENFIOOFFSET - 4)
#define dsINPUTEXISTSWSEV dsERRSEVWARNING
#define dsINPUTEXISTSWSTDMSG "WARNING: Input file '%s' does not exist, continuing without it.\n"

#define dsREADFILEFERR (dsGENFIOOFFSET - 5)
#define dsREADFILEFSEV dsERRSEVFATAL
#define dsREADFILEFSTDMSG "ERROR: Could not read file '%s'.\n"

#define dsREADFILEWERR (dsGENFIOOFFSET - 6)
#define dsREADFILEWSEV dsERRSEVWARNING
#define dsREADFILEWSTDMSG "WARNING: Could not read file '%s', continuing without it.\n"

#define dsCLOBBERFILEERR (dsGENFIOOFFSET - 7)
#define dsCLOBBERFILESEV dsERRSEVFATAL
#define dsCLOBBERFILESTDMSG "ERROR: Could not clobber file '%s'.\n"

#define dsCREATEFILEERR (dsGENFIOOFFSET - 8)
#define dsCREATEFILESEV  dsERRSEVFATAL
#define dsCREATEFILESTDMSG "ERROR: Could not create file '%s'.\n"

#define dsOUTPUTEXISTSERR (dsGENFIOOFFSET - 9)
#define dsOUTPUTEXISTSSEV dsERRSEVFATAL
#define dsOUTPUTEXISTSSTDMSG "ERROR: Could not create file '%s', it exists and clobber = no.\n"

#define dsGETSTATUSERR (dsGENFIOOFFSET - 10)
#define dsGETSTATUSSEV dsERRSEVFATAL
#define dsGETSTATUSSTDMSG "ERROR: Cannot get status of file '%s'.\n"

#define dsTMPFILEERR (dsGENFIOOFFSET - 11)
#define dsTMPFILESEV dsERRSEVFATAL
#define dsTMPFILESTDMSG "ERROR: Could not create tmp file '%s'.  Set environment variable '%s'.\n"

#define dsNOROWSERR (dsGENFIOOFFSET - 12)
#define dsNOROWSSEV dsERRSEVFATAL
#define dsNOROWSSTDMSG "ERROR: File '%s' has no rows in it.\n"


/* Data Model related errors */

#define dsDMGENERICERR (dsGENDMOFFSET - 1)
#define dsDMGENERICSEV dsERRSEVFATAL
#define dsDMGENERICSTDMSG "ERROR: An unspecified DM error has occured.\n"

#define dsDMBLOCKCREATEERR (dsGENDMOFFSET - 2)
#define dsDMBLOCKCREATESEV dsERRSEVFATAL
#define dsDMBLOCKCREATESTDMSG "ERROR: Failed to create a DM data block in dataset '%s'.  DM Error: '%s'.\n"

#define dsDMBLOCKOPENERR (dsGENDMOFFSET - 3)
#define dsDMBLOCKOPENSEV dsERRSEVFATAL
#define dsDMBLOCKOPENSTDMSG "ERROR: Failed to open the DM data block in dataset '%s'.  DM Error: '%s'.\n"

#define dsDMTABLEOPENERR (dsGENDMOFFSET - 4)
#define dsDMTABLEOPENSEV dsERRSEVFATAL
#define dsDMTABLEOPENSTDMSG "ERROR: Failed to open a DM table in dataset '%s'.  DM Error: '%s'.\n"

#define dsDMSETVALUEERR (dsGENDMOFFSET - 5)
#define dsDMSETVALUESEV dsERRSEVFATAL
#define dsDMSETVALUESTDMSG "ERROR: Failed to set the value of a DM data descriptor in dataset '%s'.  DM Error: '%s'.\n"

#define dsDMGETVALUEERR (dsGENDMOFFSET - 6)
#define dsDMGETVALUESEV dsERRSEVFATAL
#define dsDMGETVALUESTDMSG "ERROR: Failed to get the value of a DM data descriptor in dataset '%s'.  DM Error: '%s'.\n"

#define dsDMOPENIMAGEERR (dsGENDMOFFSET - 7)
#define dsDMOPENIMAGESEV dsERRSEVFATAL
#define dsDMOPENIMAGESTDMSG "ERROR: Failed to open a DM image in dataset '%s'. DM Error: '%s'.\n"
 
#define dsDMCREATEIMAGEERR (dsGENDMOFFSET - 8)
#define dsDMCREATEIMAGESEV dsERRSEVFATAL
#define dsDMCREATEIMAGESTDMSG "ERROR: Failed to create a DM image in dataset '%s'. DM Error: '%s'.\n"

#define dsDMREADIMAGEERR (dsGENDMOFFSET - 9)
#define dsDMREADIMAGESEV dsERRSEVFATAL
#define dsDMREADIMAGESTDMSG "ERROR: Failed to read a DM image in dataset '%s'. DM Error: '%s'.\n"
 
#define dsDMWRITEIMAGEERR (dsGENDMOFFSET - 10)
#define dsDMWRITEIMAGESEV dsERRSEVFATAL
#define dsDMWRITEIMAGESTDMSG "ERROR: Failed to write a DM image in dataset '%s'. DM Error: '%s'.\n"

#define dsDMCLOSEIMAGEERR (dsGENDMOFFSET - 11)
#define dsDMCLOSEIMAGESEV dsERRSEVFATAL
#define dsDMCLOSEIMAGESTDMSG "ERROR: Failed to close a DM image in dataset '%s'. DM Error: '%s'.\n"

#define dsDMSUBSPCREATEERR (dsGENDMOFFSET - 12)
#define dsDMSUBSPCREATESEV dsERRSEVFATAL
#define dsDMSUBSPCREATESTDMSG "ERROR: Could not create subspace '%s' in dataset '%s'.  DM Error: '%s'\n"

#define dsDMSSCREATECPTERR (dsGENDMOFFSET - 13)
#define dsDMSSCREATECPTSEV dsERRSEVFATAL
#define dsDMSSCREATECPTSTDMSG "ERROR: Could not create new subspace component in dataset '%s'.  DM Error: '%s'\n"

#define dsDMSSSETELEMENTERR (dsGENDMOFFSET - 14)
#define dsDMSSSETELEMENTSEV dsERRSEVFATAL
#define dsDMSSSETELEMENTSTDMSG "ERROR: Could not set element '%s' in new subspace component in dataset '%s'.  DM Error: '%s'\n"

#define dsDMSETOPENERR (dsGENDMOFFSET - 15)
#define dsDMSETOPENSEV dsERRSEVFATAL
#define dsDMSETOPENSTDMSG "ERROR: failed to open DM dataset: '%s'\n"

#define dsDMSETKERNELERR (dsGENDMOFFSET - 16)
#define dsDMSETKERNELSEV dsERRSEVFATAL
#define dsDMSETKERNELSTDMSG "ERROR:Failed to set the create kernel: '%s'\n"

#define dsDMSETCREATEERR (dsGENDMOFFSET - 17)
#define dsDMSETCREATESEV dsERRSEVFATAL
#define dsDMSETCREATESTDMSG "ERROR: failed to create DM dataset: '%s'  DM error: '%s'\n"

#define dsDMSUBSPREADERR (dsGENDMOFFSET - 18)
#define dsDMSUBSPREADSEV dsERRSEVFATAL
#define dsDMSUBSPREADSTDMSG "ERROR: Failed to read subspace '%s' in dataset '%s'.  DM error: '%s'\n"

#define dsDMSUBSPOPENERR (dsGENDMOFFSET - 19)
#define dsDMSUBSPOPENSEV dsERRSEVFATAL
#define dsDMSUBSPOPENSTDMSG "ERROR: Could not open subspace '%s' in dataset '%s'. DM Error: '%s'\n"

#define dsDMDESCRIPOPENERR (dsGENDMOFFSET - 20)
#define dsDMDESCRIPOPENSEV dsERRSEVFATAL
#define dsDMDESCRIPOPENSTDMSG "ERROR: Failed to open descriptor '%s' in dataset '%s'.  DM error: '%s'\n"

#define dsDMPUTROWERR (dsGENDMOFFSET - 21)
#define dsDMPUTROWSEV dsERRSEVFATAL
#define dsDMPUTROWSTDMSG "ERROR: Problem adding row to ouput file: '%s'.  DM Error: '%s'\n"

#define dsDMREADCOLERR (dsGENDMOFFSET - 22)
#define dsDMREADCOLSEV dsERRSEVFATAL
#define dsDMREADCOLSTDMSG "ERROR: Could not read column values for '%s'.  DM Error: '%s'\n"


/* Parameter file errors */
#define dsOPENPARAMFERR (dsGENPARAMOFFSET - 1)
#define dsOPENPARAMFSEV dsERRSEVFATAL
#define dsOPENPARAMFSTDMSG "ERROR: Could not open parameter file '%s'.\n"

#define dsOPENPARAMWERR (dsGENPARAMOFFSET - 2)
#define dsOPENPARAMWSEV dsERRSEVWARNING
#define dsOPENPARAMWSTDMSG "WARNING: Could not open parameter file '%s', continuing without it.\n"

#define dsFINDPARAMFERR (dsGENPARAMOFFSET - 3)
#define dsFINDPARAMFSEV dsERRSEVFATAL
#define dsFINDPARAMFSTDMSG "ERROR: Could not find parameter '%s' in parameter file '%s'.\n"

#define dsFINDPARAMWERR (dsGENPARAMOFFSET - 4)
#define dsFINDPARAMWSEV dsERRSEVWARNING
#define dsFINDPARAMWSTDMSG "WARNING: Could not find parameter '%s' in parameter file '%s', continuing without it.\n"

#define dsPARAMVALUEFERR (dsGENPARAMOFFSET - 5)
#define dsPARAMVALUEFSEV dsERRSEVFATAL
#define dsPARAMVALUEFSTDMSG "ERROR: Parameter '%s' has an incorrect value.\n"

#define dsPARAMVALUEWERR (dsGENPARAMOFFSET - 6)
#define dsPARAMVALUEWSEV dsERRSEVWARNING
#define dsPARAMVALUEWSTDMSG "WARNING: Parameter '%s' has an incorrect value, continuing with default value.\n"

#define dsPARAMFORMATERR (dsGENPARAMOFFSET - 7)
#define dsPARAMFORMATSEV dsERRSEVFATAL
#define dsPARAMFORMATSTDMSG "ERROR: Parameter list '%s' has bad format.\n"

#define dsPARAMNULLERR (dsGENPARAMOFFSET - 8)
#define dsPARAMNULLSEV dsERRSEVFATAL
#define dsPARAMNULLSTDMSG "ERROR: File parameter '%s' is NULL.\n"

#define dsPARAMTOOLERR (dsGENPARAMOFFSET - 9)
#define dsPARAMTOOLSEV dsERRSEVFATAL
#define dsPARAMTOOLSTDMSG "ERROR: Parameter list '%s' is too long.\n"

#define dsPARAMMISSMERR (dsGENPARAMOFFSET - 10)
#define dsPARAMMISSMSEV dsERRSEVFATAL
#define dsPARAMMISSMSTDMSG "ERROR: Mismatch in parameters '%s' and '%s'.\n"


/* Errors in using the asc fitting engine */

#define dsAFFITERR (dsGENASCFITOFFSET - 1)
#define dsAFFITSEV dsERRSEVWARNING
#define dsAFFITSTDMSG "WARNING: ascfitFit failed.\n"

#define dsAFMODELVALERR (dsGENASCFITOFFSET - 2)
#define dsAFMODELVALSEV dsERRSEVWARNING
#define dsAFMODELVALSTDMSG "WARNING: ascfitGetModelValues failed.\n"

#define dsAFPARAMERR (dsGENASCFITOFFSET - 3)
#define dsAFPARAMSEV dsERRSEVWARNING
#define dsAFPARAMSTDMSG "WARNING: Problem getting ascfit parameters.\n"

#define dsAFSETCONFERR (dsGENASCFITOFFSET - 4)
#define dsAFSETCONFSEV dsERRSEVWARNING
#define dsAFSETCONFSTDMSG "WARNING: ascfitSetConfidenceDelta failed.\n"

#define dsAFSETDATA1DERR (dsGENASCFITOFFSET - 5)
#define dsAFSETDATA1DSEV dsERRSEVWARNING
#define dsAFSETDATA1DSTDMSG "WARNING: ascfitSetData1D failed.\n"

#define dsAFSETMETHODERR (dsGENASCFITOFFSET - 6)
#define dsAFSETMETHODSEV dsERRSEVWARNING
#define dsAFSETMETHODSTDMSG "WARNING: ascfitSetMethod failed.\n"

#define dsAFSETMODELERR (dsGENASCFITOFFSET - 7)
#define dsAFSETMODELSEV dsERRSEVWARNING
#define dsAFSETMODELSTDMSG  "WARNING: ascfitSetModel failed.\n"

#define dsAFSETPARAMERR (dsGENASCFITOFFSET - 8) 
#define dsAFSETPARAMSEV dsERRSEVWARNING
#define dsAFSETPARAMSTDMSG "WARNING: ascfitSetParams failed.\n"

#define dsAFSETRANGEERR (dsGENASCFITOFFSET - 9)
#define dsAFSETRANGESEV dsERRSEVWARNING
#define dsAFSETRANGESTDMSG  "WARNING: ascfitSetRange failed.\n"

#define dsAFSETSTATERR (dsGENASCFITOFFSET - 10)
#define dsAFSETSTATSEV dsERRSEVWARNING
#define dsAFSETSTATSTDMSG "WARNING: ascfitSetStatistic failed.\n"

/* Timing errors */

#define dsGETGTIERR (dsGENTIMEOFFSET - 1)
#define dsGETGTISEV dsERRSEVFATAL
#define dsGETGTISTDMSG "ERROR: Failed to get good time interval.\n"

#define dsGETTIMEERR (dsGENTIMEOFFSET - 2)
#define dsGETTIMESEV dsERRSEVFATAL
#define dsGETTIMESTDMSG "ERROR: Failed to get start time and stop time.\n"

#define dsTIMESORTERR (dsGENTIMEOFFSET - 3)
#define dsTIMESORTSEV dsERRSEVFATAL
#define dsTIMESORTSTDMSG "ERROR: Input file '%s' is not time-sorted.\n"

#define dsTIMESORTWERR (dsGENTIMEOFFSET - 4)
#define dsTIMESORTWSEV dsERRSEVWARNING
#define dsTIMESORTWSTDMSG "WARNING: Time smaller than the previous in the following rows: '%s' \n"

#define dsFILEORDERERR (dsGENTIMEOFFSET - 5)
#define dsFILEORDERSEV dsERRSEVFATAL
#define dsFILEORDERSTDMSG "ERROR: TSTART of the current file '%s' is smaller than TSTOP of the previous file '%s' \n"  

#define dsTIMERNGERR (dsGENTIMEOFFSET - 6)
#define dsTIMERNGSEV dsERRSEVFATAL
#define dsTIMERNGSTDMSG "ERROR: Start time is after the stop time.\n"


/* stack lib errors */
#define dsSTKGENERICERR (dsGENSTACKOFFSET - 1)
#define dsSTKGENERICSEV dsERRSEVFATAL
#define dsSTKGENERICSTDMSG "ERROR: An unspecified stack library error has occurred.\n"

#define dsSTKBLDERR (dsGENSTACKOFFSET - 2)
#define dsSTKBLDSEV dsERRSEVFATAL
#define dsSTKBLDSTDMSG "ERROR: Cannot build stack from '%s'.\n"

#define dsSTKEMPTYERR (dsGENSTACKOFFSET - 3)
#define dsSTKEMPTYSEV dsERRSEVFATAL
#define dsSTKEMPTYSTDMSG "ERROR: No more files in the stack.\n"

#define dsSTKREADERR (dsGENSTACKOFFSET - 4)
#define dsSTKREADSEV dsERRSEVFATAL
#define dsSTKREADSTDMSG "ERROR: Cannot read filename #%d from stack.\n"

#define dsSTKCLOSEERR (dsGENSTACKOFFSET - 5)
#define dsSTKCLOSESEV dsERRSEVFATAL
#define dsSTKCLOSESTDMSG "ERROR: Problem closing the stack. \n"

#define dsSTKADDERR (dsGENSTACKOFFSET - 6)
#define dsSTKADDSEV dsERRSEVFATAL
#define dsSTKADDSTDMSG "ERROR: Could not add item to stack. \n"

#define dsSTKDELERR (dsGENSTACKOFFSET - 7)
#define dsSTKDELSEV dsERRSEVFATAL
#define dsSTKDELSTDMSG "ERROR: Could not delete item from stack. \n"

/* caldb lib errs */

#define dsCALDBGENFERR ( dsGENCALDBOFFSET - 1)
#define dsCALDBGENFSEV  dsERRSEVFATAL
#define dsCALDBGENFMSG "ERROR: Problems with CALDB interface.\n"

#define dsCALDBGENWERR ( dsGENCALDBOFFSET - 2)
#define dsCALDBGENWSEV  dsERRSEVWARNING
#define dsCALDBGENWMSG "WARNING: Problems with CALDB interface.\n"



