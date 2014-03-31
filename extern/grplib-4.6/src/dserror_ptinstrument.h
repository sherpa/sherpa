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
#define _DSERROR_PTINSTRUMENT_H

#define dsPTINSTRUMENTERROFFSET     -3000
#define dsPTINSTRUMENTNUMERRORS     158


/* errors for ACIS_PROCESS_EVENTS */
#define dsAPEEMPTYSTACKERR (dsPTINSTRUMENTERROFFSET - 0)
#define dsAPEEMPTYSTACKSEV dsERRSEVFATAL
#define dsAPEEMPTYSTACKSTDMSG "ERROR: The input event stack does not contain any files."

#define dsAPEGRADEFILEERR (dsPTINSTRUMENTERROFFSET - 1)
#define dsAPEGRADEFILESEV dsERRSEVFATAL
#define dsAPEGRADEFILESTDMSG "ERROR: The input grade scheme file %s could not be opened."

#define dsAPEGRADEFILESZERR (dsPTINSTRUMENTERROFFSET - 2)
#define dsAPEGRADEFILESZSEV dsERRSEVFATAL
#define dsAPEGRADEFILESZSTDMSG "ERROR: The grade scheme file %s is missing entries or contains unrecognized columns."

#define dsAPEFOAFPDNEERR (dsPTINSTRUMENTERROFFSET - 3)
#define dsAPEFOAFPDNESEV dsERRSEVWARNING
#define dsAPEFOAFPDNESTDMSG "WARNING: foafp_x, foafp_y, or foafp_z data does not exist in file %s ot performing 'winsize' corrections."
 
#define dsAPEWINSIZEERR (dsPTINSTRUMENTERROFFSET - 4)
#define dsAPEWINSIZESEV dsERRSEVFATAL
#define dsAPEWINSIZESTDMSG "ERROR: All files in the stack must have the WINSIZE settings."

#define dsAPEGAINFILEERR (dsPTINSTRUMENTERROFFSET - 5)
#define dsAPEGAINFILESEV dsERRSEVFATAL
#define dsAPEGAINFILESTDMSG "ERROR: The gain file %s does contain all of the necessary columns."

#define dsAPEDOPINOGAINERR  (dsPTINSTRUMENTERROFFSET - 6)
#define dsAPEDOPINOGAINSEV dsERRSEVFATAL
#define dsAPEDOPINOGAINSTDMSG "ERROR: The calculate_pi option is set to yes but a gain file was not provided." 

#define dsAPEBADGAINROWERR (dsPTINSTRUMENTERROFFSET - 7)
#define dsAPEBADGAINROWSEV dsERRSEVFATAL
#define dsAPEBADGAINROWSTDMSG "ERROR: A row in the gain table contains an invalid ccd id or ccd node value."

#define dsAPEGAINMISSERR (dsPTINSTRUMENTERROFFSET - 8)
#define dsAPEGAINMISSSEV dsERRSEVFATAL
#define dsAPEGAINMISSSTDMSG "ERROR: The gain file is missing one or more entries."

#define dsAPEGAINCCDRNGERR (dsPTINSTRUMENTERROFFSET - 9)
#define dsAPEGAINCCDRNGSEV dsERRSEVFATAL
#define dsAPEGAINCCDRNGSTDMSG "ERROR: An event with an out of range ccd id (%hd) was detected."

#define dsAPEGAINNODERNGERR (dsPTINSTRUMENTERROFFSET - 10)
#define dsAPEGAINNODERNGSEV dsERRSEVFATAL
#define dsAPEGAINNODERNGSTDMSG "ERROR: An event with an out of range ccd id (%hd) was detected."

#define dsAPETSTARTERR (dsPTINSTRUMENTERROFFSET - 11)
#define dsAPETSTARTSEV dsERRSEVWARNING 
#define dsAPETSTARTSTDMSG "WARNING: The input event times start at %f but the TSTART specified in %s is %f."

#define dsAPETSTOPERR (dsPTINSTRUMENTERROFFSET - 12)
#define dsAPETSTOPSEV dsERRSEVWARNING 
#define dsAPETSTOPSTDMSG "WARNING: The input event times end at %f but the TSTOP specified in %s is %f."

#define dsAPEOBSTIMERNGERR (dsPTINSTRUMENTERROFFSET - 13)
#define dsAPEOBSTIMERNGSEV dsERRSEVWARNING
#define dsAPEOBSTIMERNGSTDMSG "WARNING: The tstop time (%f) is earlier than the tstart time (%f) in %s."

#define dsAPEDMDATATYPEERR (dsPTINSTRUMENTERROFFSET - 14)
#define dsAPEDMDATATYPESEV dsERRSEVFATAL
#define dsAPEDMDATATYPESTDMSG "ERROR: The specified output column %s is composed of un unknown or unsupported data type."

#define dsAPEUNKNOWNSYSERR (dsPTINSTRUMENTERROFFSET - 15)
#define dsAPEUNKNOWNSYSSEV dsERRSEVWARNING
#define dsAPEUNKNOWNSYSSTDMSG "WARNING: Unable to set the TLMIN/TLMAX values for output event file column %s due to unknown instrument configuration."

/*Errors added 28 February 2007 for ACIS_PROCESS_EVENTS */
#define dsAPEPULSEHEIGHTERR (dsPTINSTRUMENTERROFFSET - 16)
#define dsAPEPULSEHEIGHTSEV dsERRSEVWARNING
#define dsAPEPULSEHEIGHTSTDMSG "WARNING: pulse height is less than split threshold when performing serial CTI adjustment.\n"

#define dsAPEREADERR (dsPTINSTRUMENTERROFFSET - 17)
#define dsAPEREADSEV dsERRSEVFATAL
#define dsAPEREADSTDMSG "ERROR: Problem reading %s.\n"

#define dsAPEDXYERR (dsPTINSTRUMENTERROFFSET - 18)
#define dsAPEDXYSEV dsERRSEVWARNING
#define dsAPEDXYSTDMSG "WARNING: Problem computing previous diffx/y values in CTI adjustment.\n"

#define dsAPEDIVZEROERR (dsPTINSTRUMENTERROFFSET - 19)
#define dsAPEDIVZEROSEV dsERRSEVWARNING
#define dsAPEDIVZEROSTDMSG "WARNING: Attempted to divide by zero- please verify that the gain table entries are correct.\n"

#define dsAPETARGERR (dsPTINSTRUMENTERROFFSET - 20)
#define dsAPETARGSEV dsERRSEVWARNING
#define dsAPETARGSTDMSG "WARNING: the RA_TARG and DEC_TARG coordinate fall off of the CCD during the observation.  Times may be bad.\n"

#define dsAPENUMCOLSERR (dsPTINSTRUMENTERROFFSET - 21)
#define dsAPENUMCOLSSEV dsERRSEVWARNING
#define dsAPENUMCOLSSTDMSG "WARNING: %s does not contain all required columns.\n"

#define dsAPEREADCTIFILEERR (dsPTINSTRUMENTERROFFSET - 22)
#define dsAPEREADCTIFILESEV dsERRSEVWARNING
#define dsAPeREADCTIFILESTDMSG "WARNING: problem reading ctifile, cti adjustment will not be applied.\n"

#define dsAPEREADTGAINFILEERR (dsPTINSTRUMENTERROFFSET - 23)
#define dsAPEREADTGAINFILESEV dsERRSEVWARNING
#define dsAPEREADTGAINFILESTDMSG "WARNING: problem reading tgainfile, tgain adjustment will not be applied.\n"

#define dsAPEPHAZEROERR (dsPTINSTRUMENTERROFFSET - 24)
#define dsAPEPHAZEROSEV dsERRSEVWARNING
#define dsAPEPHAZEROSTDMSG "WARNING: an event's summed pha is zero or null.\n"

#define dsAPESETTRAPDENSERR (dsPTINSTRUMENTERROFFSET - 25)
#define dsAPESETTRAPDENSSEV dsERRSEVWARNING
#define dsAPESETTRAPDENSSTDMSG "WARNING: Error setting TRAPDENS value.\n"

#define dsAPEMEMERR (dsPTINSTRUMENTERROFFSET - 26)
#define dsAPEMEMSEV dsERRSEVWARNING
#define dsAPEMEMSTDMSG "WARNING: Out of memory.\n"

#define dsAPEREADIMERR (dsPTINSTRUMENTERROFFSET - 27)
#define dsAPEREADIMSEV dsERRSEVWARNING
#define dsAPEREADIMSTDMSG "WARNING: Error reading image data.\n"

#define dsAPEIMTYPEERR (dsPTINSTRUMENTERROFFSET - 28)
#define dsAPEIMTYPESEV dsERRSEVWARNING
#define dsAPEIMTYPESTDMSG "WARNING: Unrecognized image data type.\n"

#define dsAPEIMDIMERR (dsPTINSTRUMENTERROFFSET - 29)
#define dsAPEIMDIMSEV dsERRSEVWARNING
#define dsAPEIMDIMSTDMSG "WARNING: Couldn't get image dimensions.\n"

#define dsAPEMALLOCERR (dsPTINSTRUMENTERROFFSET - 30)
#define dsAPEMALLOCSEV dsERRSEVWARNING
#define dsAPEMALLOCSTDMSG "WARNING: Malloc failed for line buffer: Out of memory.\n" 

#define dsAPEREADTRAPDENSERR (dsPTINSTRUMENTERROFFSET - 31)
#define dsAPEREADTRAPDENSSEV dsERRSEVWARNING
#define dsAPEREADTRAPDENSSTDMSG "WARNING: problem reading from CTI trapdens map.\n"

#define dsAPECCDNODEERR (dsPTINSTRUMENTERROFFSET - 32)
#define dsAPECCDNODESEV dsERRSEVWARNING
#define dsAPECCDNODESTDMSG "WARNING: An invalid ccdnode (valid range is 0-3) was converted to 0 during processing.\n"
/*END 28 Feb 2007 Update  */
 


/* errors for ACIS_FORMAT_EVENTS */
#define dsAFESTACKSIZEERR (dsPTINSTRUMENTERROFFSET - 50)
#define dsAFESTACKSIZESEV dsERRSEVFATAL
#define dsAFESTACKSIZESTDMSG "ERROR: The input bias and exposure stacks must either be empty (NONE) or contain the same number of elements as the input event stack."
 
#define dsAFEEXPSTATSERR (dsPTINSTRUMENTERROFFSET - 51)
#define dsAFEEXPSTATSSEV dsERRSEVFATAL
#define dsAFEEXPSTATSSTDMSG "ERROR: Can not create exposure statistics file %s because the input exposure file stack is empty."
 
#define dsAFEBIASLOADERR (dsPTINSTRUMENTERROFFSET - 52)
#define dsAFEBIASLOADSEV dsERRSEVFATAL
#define dsAFEBIASLOADSTDMSG "ERROR: Unable to load bias map %s."
 
#define dsAFENOBIASDATAERR (dsPTINSTRUMENTERROFFSET - 53)
#define dsAFENOBIASDATASEV dsERRSEVFATAL
#define dsAFENOBIASDATASTDMSG "ERROR: A bias map must be provided for all modes except FAINT w/ BIAS"
 
#define dsAFEGRDMDBIASERR (dsPTINSTRUMENTERROFFSET - 54)
#define dsAFEGRDMDBIASSEV dsERRSEVFATAL
#define dsAFEGRDMDBIASSTDMSG "ERROR: Unable to perform bias corrections on graded mode data."
 
#define dsAFEUNKNOWNMODEERR (dsPTINSTRUMENTERROFFSET - 55)
#define dsAFEUNKNOWNMODESEV dsERRSEVFATAL
#define dsAFEUNKNOWNMODESTDMSG "ERROR: Unable to determine the data type of %s from the READMODE and DATAMODE keywords."
 
#define dsAFEMIXEDMODEERR (dsPTINSTRUMENTERROFFSET - 56)
#define dsAFEMIXEDMODESEV dsERRSEVFATAL
#define dsAFEMIXEDMODESTDMSG "ERROR: The readmode and datamode values of all input event files in the stack must be identical."
 
#define dsAFEBIASDIFFERERR (dsPTINSTRUMENTERROFFSET - 57)
#define dsAFEBIASDIFFERSEV dsERRSEVWARNING
#define dsAFEBIASDIFFERSTDMSG "WARNING: A bias pixel value in the output bias map has changed value between different events occurring in the same area."
 
#define dsAFENOBIASCORRERR (dsPTINSTRUMENTERROFFSET - 58)
#define dsAFENOBIASCORRSEV dsERRSEVWARNING
#define dsAFENOBIASCORRSTDMSG "WARNING: Not performing bias correction of faint mode data because bias_correct parameter is set to no."
 
#define dsAFEOCCHIPRANGEERR (dsPTINSTRUMENTERROFFSET - 59)
#define dsAFEOCCHIPRANGESEV dsERRSEVWARNING
#define dsAFEOCCHIPRANGESTDMSG "WARNING: Event chip position is outside of valid range. Overclock corrections not applied to this event."

#define dsAFEBSCHIPRANGEERR (dsPTINSTRUMENTERROFFSET - 60)
#define dsAFEBSCHIPRANGESEV dsERRSEVWARNING
#define dsAFEBSCHIPRANGESTDMSG "WARNING: Event chip position is outside of valid range. Bias values for this event were not applied to the output bias map."

#define dsAFECCDRANGEERR (dsPTINSTRUMENTERROFFSET - 61)
#define dsAFECCDRANGESEV dsERRSEVFATAL
#define dsAFECCDRANGESTDMSG "ERROR: The default ccd value of %hd obtained from the input file header is out of valid range (0-9)." 

#define dsAFEBADPMODEERR (dsPTINSTRUMENTERROFFSET - 62)
#define dsAFEBADPMODESEV dsERRSEVWARNING
#define dsAFEBADPMODESTDMSG "WARNING: The specified acis mode is not recognized- can not perform bad pixel check." 

#define dsAFEBPCHIPRANGEERR (dsPTINSTRUMENTERROFFSET - 63)
#define dsAFEBPCHIPRANGESEV dsERRSEVWARNING
#define dsAFEBPCHIPRANGESTDMSG "WARNING: Event chip position is outside of valid range. Bad pixel checking was not performed for this event."
 
#define dsAFEBADPCNTERR (dsPTINSTRUMENTERROFFSET - 64)
#define dsAFEBADPCNTSEV dsERRSEVWARNING
#define dsAFEBADPCNTSTDMSG "WARNING: Event island contains 1 or more bad pixels."

#define dsAFEBADPPOSERR (dsPTINSTRUMENTERROFFSET - 65)
#define dsAFEBADPPOSSEV dsERRSEVWARNING
#define dsAFEBADPPOSSTDMSG "WARNING: A specified bad pixel/column was not used in the bad pixel check because it contains an invalid point."

/* ACIS Collate Events Errors */
#define dsACEDTYCYCLE0ERR (dsPTINSTRUMENTERROFFSET - 200)
#define dsACEDTYCYCLE0SEV dsERRSEVFATAL
#define dsACEDTYCYCLE0STDMSG "ERROR: Dtycycle of parameter block is 0. Event and exposure files are not in interleaved mode."
 
#define dsACEDTYCYCRNGERR (dsPTINSTRUMENTERROFFSET - 201)
#define dsACEDTYCYCRNGSEV dsERRSEVWARNING
#define dsACEDTYCYCRNGSTDMSG "WARNING: Dtycycle value of %d is out of valid range (0..15)."
 
#define dsACEEVTEXPSTKSZERR (dsPTINSTRUMENTERROFFSET - 202)
#define dsACEEVTEXPSTKSZSEV dsERRSEVFATAL
#define dsACEEVTEXPSTKSZSTDMSG "ERROR: If exposure and event stacks both contain elements then they must contain the same number of elements."
 
#define dsACEINOUTSTKSZERR (dsPTINSTRUMENTERROFFSET - 203)
#define dsACEINOUTSTKSZSEV dsERRSEVFATAL
#define dsACEINOUTSTKSZSTDMSG "ERROR: The output stacks must contain the same number of elements as the input stacks."
 
#define dsACECOLUMNAMEERR (dsPTINSTRUMENTERROFFSET - 204)
#define dsACECOLUMNAMESEV dsERRSEVFATAL
#define dsACECOLUMNAMESTDMSG "ERROR: One or more columns in file %s not recognized as a valid level 1 column name."
 
#define dsACEEXPCOLUMNAMEERR (dsPTINSTRUMENTERROFFSET - 205)
#define dsACEEXPCOLUMNAMESEV dsERRSEVFATAL
#define dsACEEXPCOLUMNAMESTDMSG "ERROR: Column %s in file %s not recognized as a valid level 1 exposure file column name."
 
#define dsACEWRITEEVTERR (dsPTINSTRUMENTERROFFSET - 206)
#define dsACEWRITEEVTSEV dsERRSEVFATAL
#define dsACEWRITEEVTSTDMSG "ERROR: Specified output event file column either contains an unrecognized column name or column datatype."
 
#define dsACEWRITEEXPERR (dsPTINSTRUMENTERROFFSET - 207)
#define dsACEWRITEEXPSEV dsERRSEVFATAL
#define dsACEWRITEEXPSTDMSG "ERROR: Specified output column in exposure file %s either contains an unrecognized column name or column datatype."
 
#define dsACEDTYCYCMISSERR (dsPTINSTRUMENTERROFFSET - 208)
#define dsACEDTYCYCMISSSEV dsERRSEVFATAL
#define dsACEDTYCYCMISSSTDMSG "ERROR: Dtycycle keyord is missing from the parameter block file."
 
#define dsACEEXPTIMEMISSERR (dsPTINSTRUMENTERROFFSET - 209)
#define dsACEEXPTIMEMISSSEV dsERRSEVWARNING
#define dsACEEXPTIMEMISSSTDMSG "WARNING: %s keyord is missing from the parameter block file."
 
#define dsACEEXPTIMERNGERR (dsPTINSTRUMENTERROFFSET - 210)
#define dsACEEXPTIMERNGSEV dsERRSEVWARNING
#define dsACEEXPTIMERNGSTDMSG "WARNING: %s value of %f is out of valid range (0.0 - 100.0 in tenths of seconds)."

#define dsACEPBKEXTRAERR (dsPTINSTRUMENTERROFFSET - 211)
#define dsACEPBKEXTRASEV dsERRSEVWARNING
#define dsACEPBKEXTRASTDMSG "WARNING: More than 1 parameter block supplied. Only using the first one."

#define dsACEPBKNONEERR (dsPTINSTRUMENTERROFFSET - 212)
#define dsACEPBKNONESEV dsERRSEVFATAL
#define dsACEPBKNONESTDMSG "Error: No parameter blocks supplied in file %s."


/* errors for ACIS_BUILD_CHIP_GTI */
#define dsACISGTINAMESERR (dsPTINSTRUMENTERROFFSET - 250)
#define dsACISGTINAMESSEV dsERRSEVFATAL
#define dsACISGTINAMESSTDMSG "ERROR: The value of keyword %s, %s, in file %s is not an allowed value.  See help file.\n"

#define dsACISGTIDTCERR (dsPTINSTRUMENTERROFFSET - 251)
#define dsACISGTIDTCSEV dsERRSEVWARNING
#define dsACISGTIDTCSTDMSG "WARNING: Calculated Dead Time Correction was %g, setting it to 1.0.\n"

#define dsACISGTIROWERR (dsPTINSTRUMENTERROFFSET - 252)
#define dsACISGTIROWSEV dsERRSEVWARNING
#define dsACISGTIROWSTDMSG "WARNING: The %s extension in file %s does not have enough rows to determine the record length.\nLooking instead for %s keyword to get record length from.\n"


/* ACIS_MERGE_GTI ERRORS */
#define dsAMGNOBLOCKERR (dsPTINSTRUMENTERROFFSET - 300)
#define dsAMGNOBLOCKSEV dsERRSEVWARNING
#define dsAMGNOBLOCKSTDMSG "WARNING: No Block specified, attempting to open block %s in file %s.\n"



/* DESTREAK errors (ACIS filtering tool) */
#define dsDSINVALIDMAXERR (dsPTINSTRUMENTERROFFSET - 325)
#define dsDSINVALIDMAXSEV dsERRSEVFATAL
#define dsDSINVALIDMAXSTDMSG "ERROR: Invalid max:  %s\n"

#define dsDSEVENTQEMPTYERR (dsPTINSTRUMENTERROFFSET - 326)
#define dsDSEVENTQEMPTYSEV dsERRSEVFATAL
#define dsDSEVENTQEMPTYSTDMSG "ERROR: The event queue from the input table is empty.  Unable to continue.\n"

#define dsDSNOQORMAPFUNERR (dsPTINSTRUMENTERROFFSET - 327)
#define dsDSNOQORMAPFUNSEV dsERRSEVFATAL
#define dsDSNOQORMAPFUNSTDMSG "ERROR: There is no event queue or mapped function to perform.  Unable to continue.\n"

#define dsDSNOMOREROWSERR (dsPTINSTRUMENTERROFFSET - 328)
#define dsDSNOMOREROWSSEV dsERRSEVFATAL
#define dsDSNOMOREROWSSTDMSG "INTERNAL ERROR: no more rows in file at row %ld\n"

#define dsDSFILEOPENERR (dsPTINSTRUMENTERROFFSET - 329)
#define dsDSFILEOPENSEV dsERRSEVFATAL
#define dsDSFILEOPENSTDMSG "ERROR:  Failed opening %s for writing.\n"

#define dsDSNOINPUTERR (dsPTINSTRUMENTERROFFSET - 330)
#define dsDSNOINPUTSEV dsERRSEVFATAL
#define dsDSNOINPUTSTDMSG "ERROR:  No input filename specified.\n"

#define dsDSNOOUTPUTERR (dsPTINSTRUMENTERROFFSET - 331)
#define dsDSNOOUTPUTSEV dsERRSEVFATAL
#define dsDSNOOUTPUTSTDMSG "ERROR:  No output filename specified.\n"

#define dsDSMISSINGCOLSERR (dsPTINSTRUMENTERROFFSET - 332)
#define dsDSMISSINGCOLSSEV dsERRSEVFATAL
#define dsDSMISSINGCOLSSTDMSG "ERROR:  Event file does not have all required columns.\n"

#define dsDSDOMALLOCERR (dsPTINSTRUMENTERROFFSET - 333)
#define dsDSDOMALLOCSEV dsERRSEVFATAL
#define dsDSDOMALLOCSTDMSG "ERROR:  do_malloc failed, event not initialized.\n"

#define dsDSQREINITERR (dsPTINSTRUMENTERROFFSET - 334)
#define dsDSQREINITSEV dsERRSEVFATAL
#define dsDSQREINITSTDMSG "ERROR:  unable to reinitialize queue, bad parameters.\n"

#define dsDSDOREALLOCERR (dsPTINSTRUMENTERROFFSET - 335)
#define dsDSDOREALLOCSEV dsERRSEVFATAL
#define dsDSDOREALLOCSTDMSG "ERROR:  do_realloc failed, event queue failed to resize.\n"

#define dsDSNOFITERR (dsPTINSTRUMENTERROFFSET - 336)
#define dsDSNOFITSEV dsERRSEVFATAL
#define dsDSNOFITSTDMSG "ERROR:  Failed fitting row-count distribution for ccd %ld node %ld\n"

#define dsDSQOVERFLOWERR (dsPTINSTRUMENTERROFFSET - 337)
#define dsDSQOVERFLOWSEV dsERRSEVWARNING
#define dsDSQOVERFLOWSTDMSG "WARNING:  Queue overflowed with >%ld events in a single frame\n"

#define dsDSQNEWSIZEERR (dsPTINSTRUMENTERROFFSET - 338)
#define dsDSQNEWSIZESEV dsERRSEVWARNING
#define dsDSQNEWSIZESTDMSG "WARNING:  Reallocating with bufsize= %ld\n"

#define dsDSNOPROCFRAMEERR (dsPTINSTRUMENTERROFFSET - 339)
#define dsDSNOPROCFRAMESEV dsERRSEVWARNING
#define dsDSNOPROCFRAMESTDMSG "ERROR:  Failed processing frame %ld\n"

#define dsDSNOSKYCOLERR (dsPTINSTRUMENTERROFFSET - 340)
#define dsDSNOSKYCOLSEV dsERRSEVWARNING
#define dsDSNOSKYCOLSTDMSG "WARNING:  No sky column found, mask will be ignored\n"

#define dsDSNOCCDCOLERR (dsPTINSTRUMENTERROFFSET - 341)
#define dsDSNOCCDCOLSEV dsERRSEVWARNING
#define dsDSNOCCDCOLSTDMSG "WARNING:  No ccd_id column found, all events assigned to ccd 0\n"

#define dsDSTOOMANYCHIPSERR (dsPTINSTRUMENTERROFFSET - 342)
#define dsDSTOOMANYCHIPSSEV dsERRSEVWARNING
#define dsDSTOOMANYCHIPSSTDMSG "WARNING:  Instrument has too many chips: %ld"

#define dsDSTOOMANYNODESERR (dsPTINSTRUMENTERROFFSET - 343)
#define dsDSTOOMANYNODESSEV dsERRSEVWARNING
#define dsDSTOOMANYNODESSTDMSG "WARNING:  Instrument has too many nodes: %ld"

#define dsDSNONNEGSLOPEERR (dsPTINSTRUMENTERROFFSET - 344)
#define dsDSNONNEGSLOPESEV dsERRSEVWARNING
#define dsDSNONNEGSLOPESTDMSG "WARNING:  Non-negative slope of fit.\n  Possibly due to low number of counts.\n  CCD %i Node %i slope %f\n"

#define dsDSNONODEIDCOLERR (dsPTINSTRUMENTERROFFSET - 345)
#define dsDSNONODEIDCOLSEV dsERRSEVWARNING
#define dsDSNONODEIDCOLSTDMSG "WARNING:  No node_id column found, all events assigned to node 0.\n"

#define dsDSNOEXPTIMEERR (dsPTINSTRUMENTERROFFSET - 346)
#define dsDSNOEXPTIMESEV dsERRSEVWARNING
#define dsDSNOEXPTIMESTDMSG "WARNING:  No positive value for exptime specified.\n    EXPTIME keyword not found.  No timefile will be written.\n"

#define dsDSBADCCDIDERR (dsPTINSTRUMENTERROFFSET - 347)
#define dsDSBADCCDIDSEV dsERRSEVWARNING
#define dsDSBADCCDIDSTDMSG "WARNING:  Improper ccd_id parameter.  Filtering all chips."

#define dsDSNODETNAMERR (dsPTINSTRUMENTERROFFSET - 348)
#define dsDSNODETNAMSEV dsERRSEVWARNING
#define dsDSNODETNAMSTDMSG "WARNING: DETNAM keyword missing in infile header."



/* HRC_PROCESS_EVENTS errors */
#define dsHPEBADCOORDRNGERR (dsPTINSTRUMENTERROFFSET - 400)
#define dsHPEBADCOORDRNGSEV dsERRSEVFATAL
#define dsHPEBADCOORDRNGSTDMSG "ERROR: The coordinate transformation starting point specified in the .par file is invalid."

#define dsHPETIMEORDERERR (dsPTINSTRUMENTERROFFSET - 401)
#define dsHPETIMEORDERSEV dsERRSEVFATAL
#define dsHPETIMEORDERSTDMSG "ERROR: Event files must be in chronological order."

#define dsHPESTARTSTOPERR (dsPTINSTRUMENTERROFFSET - 402)
#define dsHPESTARTSTOPSEV dsERRSEVFATAL
#define dsHPESTARTSTOPSTDMSG "ERROR: Event file TSTOP time must be greater than TSTART time."

#define dsHPEDEGAPALLOCERR (dsPTINSTRUMENTERROFFSET - 403)
#define dsHPEDEGAPALLOCSEV dsERRSEVFATAL
#define dsHPEDEGAPALLOCSTDMSG "ERROR: Degap table memory allocation failed."

#define dsHPEDEGAPLOADERR (dsPTINSTRUMENTERROFFSET - 404)
#define dsHPEDEGAPLOADSEV dsERRSEVFATAL
#define dsHPEDEGAPLOADSTDMSG "ERROR: Unable to load degap file into memory."

#define dsHPEADCALLOCERR (dsPTINSTRUMENTERROFFSET - 405)
#define dsHPEADCALLOCSEV dsERRSEVFATAL
#define dsHPEADCALLOCSTDMSG "ERROR: ADC Correction table memory allocation failed."

#define dsHPEADCLOADERR (dsPTINSTRUMENTERROFFSET - 406)
#define dsHPEADCLOADSEV dsERRSEVFATAL
#define dsHPEADCLOADSTDMSG "ERROR: Unable to load ADC Correction file into memory."

#define dsHPEOUTCOLUMNERR (dsPTINSTRUMENTERROFFSET - 407) 
#define dsHPEOUTCOLUMNSEV dsERRSEVFATAL
#define dsHPEOUTCOLUMNSTDMSG "ERROR: Coordinate columns specified in the output eventdef do not correspond to requested coordinate transformations."  

#define dsHPEALIGNMENTERR (dsPTINSTRUMENTERROFFSET - 408)
#define dsHPEALIGNMENTSEV dsERRSEVFATAL
#define dsHPEALIGNMENTSTDMSG "ERROR: The alignment file could not be opened or is missing required columns." 

#define dsHPEASPECTERR (dsPTINSTRUMENTERROFFSET - 409)
#define dsHPEASPECTSEV dsERRSEVFATAL
#define dsHPEASPECTSTDMSG "ERROR: The aspect file could not be opened or is missing required columns." 

#define dsHPEEVENTSEQERR (dsPTINSTRUMENTERROFFSET - 410)
#define dsHPEEVENTSEQSEV dsERRSEVWARNING
#define dsHPEEVENTSEQSTDMSG "WARNING: Out of sequence events discovered in %s."

#define dsHPEALLBADFILESERR (dsPTINSTRUMENTERROFFSET - 411)
#define dsHPEALLBADFILESSEV dsERRSEVFATAL
#define dsHPEALLBADFILESSTDMSG "ERROR: Unable to successfully process any input event files."

#define dsHPEFILEEXISTSERR (dsPTINSTRUMENTERROFFSET - 412)
#define dsHPEFILEEXISTSSEV dsERRSEVFATAL
#define dsHPEFILEEXISTSSTDMSG "ERROR: Could not create output file %s because the file already exists."

#define dsHPEWRITEEVTERR (dsPTINSTRUMENTERROFFSET - 413)
#define dsHPEWRITEEVTSEV dsERRSEVWARNING
#define dsHPEWRITEEVTSTDMSG "WARNING: Could not write event data to unrecognized column of %s." 

#define dsHPEINSTRUMEPARERR  (dsPTINSTRUMENTERROFFSET - 414)
#define dsHPEINSTRUMEPARSEV dsERRSEVFATAL
#define dsHPEINSTRUMEPARSTDMSG "ERROR: Unable to open the file specified by %s parameter in %s."

#define dsHPEBADEVTFILEERR (dsPTINSTRUMENTERROFFSET - 415)
#define dsHPEBADEVTFILESEV dsERRSEVWARNING  
#define dsHPEBADEVTFILESTDMSG "WARNING: Event written to bad event file %s."  

#define dsHPELOGOPENERR (dsPTINSTRUMENTERROFFSET - 416)
#define dsHPELOGOPENSEV dsERRSEVWARNING  
#define dsHPELOGOPENSTDMSG "WARNING: Could not open logfile %s. Writing log to stdout." 

#define dsHPESTAGEDNEERR (dsPTINSTRUMENTERROFFSET - 417)
#define dsHPESTAGEDNESEV dsERRSEVWARNING
#define dsHPESTAGEDNESTDMSG "WARNING: Stage keywords (STF_X/Y/Z) not found in %s. Using zeros as the default values."

#define dsHPESTAGEANGDNEERR (dsPTINSTRUMENTERROFFSET - 418)
#define dsHPESTAGEANGDNESEV dsERRSEVWARNING
#define dsHPESTAGEANGDNESTDMSG "WARNING: Stage angle keywords (STF_ANG1/2/3) not found in %s. Using zeros as the default values."  

#define dsHPEHPYDNEERR (dsPTINSTRUMENTERROFFSET - 419)
#define dsHPEHPYDNESEV dsERRSEVWARNING
#define dsHPEHPYDNESTDMSG "WARNING: HRMA Pitch and Yaw (HRMA_PIT/YAW) keywords not found in %s. Using zeros as the default values." 

#define dsHPESETAIMPOINTERR (dsPTINSTRUMENTERROFFSET - 420)
#define dsHPESETAIMPOINTSEV dsERRSEVWARNING
#define dsHPESETAIMPOINTSTDMSG "WARNING: Sim adjustment keywords (SIM_X/Y/Z) not found in %s. Using the instrument's default aimpoint position." 

#define dsHPEDEPENDENCYERR (dsPTINSTRUMENTERROFFSET - 421)
#define dsHPEDEPENDENCYSEV dsERRSEVFATAL
#define dsHPEDEPENDENCYSTDMSG "ERROR: A data dependency on %s was not met."

#define dsHPEBADOUTCOLERR (dsPTINSTRUMENTERROFFSET - 422)
#define dsHPEBADOUTCOLSEV dsERRSEVFATAL
#define dsHPEBADOUTCOLSTDMSG "ERROR: The column %s is not recognized as a valid level 1 output event column."

#define dsHPEADCOPENERR (dsPTINSTRUMENTERROFFSET - 423) 
#define dsHPEADCOPENSEV dsERRSEVFATAL
#define dsHPEADCOPENSTDMSG "ERROR: The file %s either does not exist or is not accessible."
 
#define dsHPEADCOPENMEMERR (dsPTINSTRUMENTERROFFSET - 424) 
#define dsHPEADCOPENMEMSEV dsERRSEVFATAL
#define dsHPEADCOPENMEMSTDMSG "ERROR: Unable to allocate memory necessary to open file %s."
 
#define dsHPEADCROWCNTERR (dsPTINSTRUMENTERROFFSET - 425) 
#define dsHPEADCROWCNTSEV dsERRSEVWARNING
#define dsHPEADCROWCNTSTDMSG "WARNING: %s does not contain the expected number of rows"
 
#define dsHPEADCMISSROWERR (dsPTINSTRUMENTERROFFSET - 426) 
#define dsHPEADCMISSROWSEV dsERRSEVWARNING
#define dsHPEADCMISSROWSTDMSG "WARNING: %s does not contain values for %c tap %hd. Using 0 for P and 1 for Q."

#define dsHPEADCBADSYSERR (dsPTINSTRUMENTERROFFSET - 427) 
#define dsHPEADCBADSYSSEV dsERRSEVFATAL
#define dsHPEADCBADSYSSTDMSG "ERROR: The specified instrument %s is unrecognized or invalid."
 

/* HRC_CALC_DEAD_TIME ERRORS */
#define dsHCDTSSPSOVERLAPERR (dsPTINSTRUMENTERROFFSET - 450)
#define dsHCDTSSPSOVERLAPSEV dsERRSEVWARNING
#define dsHCDTSSPSOVERLAPSTDMSG "WARNING: The secondary and primary science do not overlap.  No records will be processed.\n"

#define dsHCDEVENNUMERR (dsPTINSTRUMENTERROFFSET - 451)
#define dsHCDEVENNUMSEV dsERRSEVWARNING
#define dsHCDEVENNUMSTDMSG "WARNING: There is not an even number of records in the nth\n\tfile in the input stack, where n is %d.  Processing will stop with the last pair of records in that file.\n"


/* HRC_BUILD_BADPIX ERRORS */
#define dsHBBDESCNOTSETERR (dsPTINSTRUMENTERROFFSET - 550)
#define dsHBBDESCNOTSETSEV dsERRSEVFATAL
#define dsHBBDESCNOTSETSTDMSG "ERROR: Unable to process the %s badpixel file because the file descriptors have not been properly set up."

#define dsHBBINVALIDSHAPEERR (dsPTINSTRUMENTERROFFSET - 551)
#define dsHBBINVALIDSHAPESEV dsERRSEVWARNING
#define dsHBBINVALIDSHAPESTDMSG "WARNING: The shape %s is currently not supported. Ignoring this element"

#define dsHBBEMPTYINFILEERR (dsPTINSTRUMENTERROFFSET - 552)
#define dsHBBEMPTYINFILESEV dsERRSEVWARNING
#define dsHBBEMPTYINFILESTDMSG "WARNING: The input bad pixel file contains no entries."

#define dsHBBOBSREREADERR (dsPTINSTRUMENTERROFFSET - 553)
#define dsHBBOBSREREADSEV dsERRSEVWARNING
#define dsHBBOBSREREADSTDMSG "WARNING: Obs.par file appears to have already been read- using previous data."

#define dsHBBNOMEMALLOCERR (dsPTINSTRUMENTERROFFSET - 554)
#define dsHBBNOMEMALLOCSEV dsERRSEVFATAL
#define dsHBBNOMEMALLOCSTDMSG "ERROR: No memory was allocated for the data structure %s."

#define dsHBBDEGAPREREADERR (dsPTINSTRUMENTERROFFSET - 555)
#define dsHBBDEGAPREREADSEV dsERRSEVWARNING
#define dsHBBDEGAPREREADSTDMSG "WARNING: The degap file appears to have already been loaded- using previous data."

#define dsHBBSETTLMINERR (dsPTINSTRUMENTERROFFSET - 556)
#define dsHBBSETTLMINSEV dsERRSEVWARNING
#define dsHBBSETTLMINSTDMSG "WARNING: chip coord. is outside of valid range. Set to TLMIN."

#define dsHBBSETTLMAXERR (dsPTINSTRUMENTERROFFSET - 557)
#define dsHBBSETTLMAXSEV dsERRSEVWARNING
#define dsHBBSETTLMAXSTDMSG "WARNING: chip coord. is outside of valid range. Set to TLMAX."


/* errors for hrc_merge_times */
#define dsHMTSSOBSTIMEERR (dsPTINSTRUMENTERROFFSET - 650)
#define dsHMTSSOBSTIMESEV dsERRSEVWARNING
#define dsHMTSSOBSTIMESTDMSG "WARNING: No secondary science records in %s overlap [TSTART=%f, TSTOP=%f] in the obspar file, %s.  No output is produced. \n"
			
#define dsHMTHRCHKOBSTIMEERR (dsPTINSTRUMENTERROFFSET - 651)
#define dsHMTHRCHKOBSTIMESEV dsERRSEVWARNING
#define dsHMTHRCHKOBSTIMESTDMSG "WARNING: No records in the hrchk_look file, %s, overlap [TSTART=%f, TSTOP=%f] in the obspar file,%s. \n"


/* errors for TG_RESOLVE_EVENTS */
#define dsTREREGLOADERR (dsPTINSTRUMENTERROFFSET - 850)
#define dsTREREGLOADSEV dsERRSEVFATAL
#define dsTREREGLOADSTDMSG "ERROR: Unable to parse or load region file %s."

#define dsTREINDATATYPEERR (dsPTINSTRUMENTERROFFSET - 851)
#define dsTREINDATATYPESEV dsERRSEVWARNING
#define dsTREINDATATYPESTDMSG "WARNING: Data value for column %s not loaded due to inability to determine input column data type."

#define dsTREOUTDATATYPEERR (dsPTINSTRUMENTERROFFSET - 852)
#define dsTREOUTDATATYPESEV dsERRSEVWARNING
#define dsTREOUTDATATYPESTDMSG "WARNING: Data value for column %s not written to output file due to inability to determine output data column type."

#define dsTRECASTFLOWERR (dsPTINSTRUMENTERROFFSET - 853)
#define dsTRECASTFLOWSEV dsERRSEVWARNING
#define dsTRECASTFLOWSTDMSG "WARNING: A potential casting overflow has occurred for data in column %s." 

#define dsTREUNKNOWNINCOLERR (dsPTINSTRUMENTERROFFSET - 854)
#define dsTREUNKNOWNINCOLSEV dsERRSEVWARNING
#define dsTREUNKNOWNINCOLSTDMSG "WARNING: Not loading data from unrecognized level 1.5 input column."

#define dsTREUNKNOWNOUTCOLERR (dsPTINSTRUMENTERROFFSET - 855)
#define dsTREUNKNOWNOUTCOLSEV dsERRSEVFATAL
#define dsTREUNKNOWNOUTCOLSTDMSG "ERROR: Could not write out unrecognized level 1.5 output column."

#define dsTREZOFPC2CHIPERR (dsPTINSTRUMENTERROFFSET - 856)
#define dsTREZOFPC2CHIPSEV dsERRSEVWARNING
#define dsTREZOFPC2CHIPSTDMSG "WARNING: Unable to determine the aspect corrected zero order focal plane position of source %hd."

#define dsTRERMENERGYFINDERR (dsPTINSTRUMENTERROFFSET - 857)
#define dsTRERMENERGYFINDSEV dsERRSEVWARNING
#define dsTRERMENERGYFINDSTDMSG "WARNING: Could not find the energy range in the rm table that matches the event."

#define dsTREBADTELESCOPERR (dsPTINSTRUMENTERROFFSET - 858)
#define dsTREBADTELESCOPSEV dsERRSEVFATAL
#define dsTREBADTELESCOPSTDMSG "ERROR: The TELESCOP keyword in %s is either invalid or differs from other files in the input event stack ."

#define dsTREMISSTELESCOPERR (dsPTINSTRUMENTERROFFSET - 859)
#define dsTREMISSTELESCOPSEV dsERRSEVWARNING
#define dsTREMISSTELESCOPSTDMSG "WARNING: The TELESCOP keyword is missing from the principal extension of file %s. Using the Axaf flight configuration as the default."

/* TG MASK SPECIFIC ERRORS */
#define dsTGMUNEQUALSTACKERR  (dsPTINSTRUMENTERROFFSET - 900)
#define dsTGMUNEQUALSTACKSEV dsERRSEVFATAL
#define dsTGMUNEQUALSTACKSTDMSG "ERROR: Input/output stacks are of unequal size.\n"

#define dsTGMINTERPOLATIONERR (dsPTINSTRUMENTERROFFSET - 901)
#define dsTGMINTERPOLATIONSEV dsERRSEVFATAL
#define dsTGMINTERPOLATIONSTDMSG "ERROR: Unable to interpolate fwhm for offaxis=%f.\n"

#define dsTGMHZEROORDERERR (dsPTINSTRUMENTERROFFSET - 902)
#define dsTGMHZEROORDERSEV dsERRSEVFATAL
#define dsTGMHZEROORDERSTDMSG "ERROR: Heg zero order radius calculations.\n"

#define dsTGMMZEROORDERERR (dsPTINSTRUMENTERROFFSET - 903)
#define dsTGMMZEROORDERSEV dsERRSEVFATAL
#define dsTGMMZEROORDERSTDMSG "ERROR: Meg zero order radius calculations.\n"

#define dsTGMLZEROORDERERR (dsPTINSTRUMENTERROFFSET - 904)
#define dsTGMLZEROORDERSEV dsERRSEVFATAL
#define dsTGMLZEROORDERSTDMSG "ERROR: Leg zero order radius calculations.\n"

#define dsTGMHWIDTHERR (dsPTINSTRUMENTERROFFSET - 905)
#define dsTGMHWIDTHSEV dsERRSEVFATAL
#define dsTGMHWIDTHSTDMSG "ERROR: Source %d is missing HEG width.\n"

#define dsTGMMWIDTHERR (dsPTINSTRUMENTERROFFSET - 906)
#define dsTGMMWIDTHSEV dsERRSEVFATAL
#define dsTGMMWIDTHSTDMSG "ERROR: Source %d is missing MEG width.\n"

#define dsTGMLWIDTHERR (dsPTINSTRUMENTERROFFSET - 907)
#define dsTGMLWIDTHSEV dsERRSEVFATAL
#define dsTGMLWIDTHSTDMSG "ERROR: Source %d is missing LEG width.\n"

#define dsTGMPIXLIBERR (dsPTINSTRUMENTERROFFSET - 908)
#define dsTGMPIXLIBSEV dsERRSEVFATAL
#define dsTGMPIXLIBSTDMSG "ERROR: PIXLIB value not found for %s\n"

#define dsTGMLENGTHCALCERR (dsPTINSTRUMENTERROFFSET - 909)
#define dsTGMLENGTHCALCSEV dsERRSEVFATAL
#define dsTGMLENGTHCALCSTDMSG "ERROR: Mask length calculation - %s\n"

#define dsTGMGRATINGMATCHERR (dsPTINSTRUMENTERROFFSET - 910)
#define dsTGMGRATINGMATCHSEV dsERRSEVFATAL
#define dsTGMGRATINGMATCHSTDMSG "ERROR: Grating name %s, not matched to known types\n"

#define dsTGMZEROYVALERR (dsPTINSTRUMENTERROFFSET - 911)
#define dsTGMZEROYVALSEV dsERRSEVFATAL
#define dsTGMZEROYVALSTDMSG "ERROR: Source %d is missing zero order y value\n"

#define dsTGMZERORADIUSERR (dsPTINSTRUMENTERROFFSET - 912)
#define dsTGMZERORADIUSSEV dsERRSEVFATAL
#define dsTGMZERORADIUSSTDMSG "ERROR: Source %d is missing zero order radius\n"

#define dsTGMHCALCWIDERR (dsPTINSTRUMENTERROFFSET - 913)
#define dsTGMHCALCWIDSEV dsERRSEVFATAL
#define dsTGMHCALCWIDSTDMSG "ERROR: Heg mask width calculations\n"

#define dsTGMMCALCWIDERR (dsPTINSTRUMENTERROFFSET - 914)
#define dsTGMMCALCWIDSEV dsERRSEVFATAL
#define dsTGMMCALCWIDSTDMSG "ERROR: Meg mask width calculations\n"

#define dsTGMLCALCWIDERR (dsPTINSTRUMENTERROFFSET - 915)
#define dsTGMLCALCWIDSEV dsERRSEVFATAL
#define dsTGMLCALCWIDSTDMSG "ERROR: Leg mask width calculations\n"

#define dsTGMUSERMASKPARSERR (dsPTINSTRUMENTERROFFSET - 916)
#define dsTGMUSERMASKPARSSEV dsERRSEVFATAL
#define dsTGMUSERMASKPARSSTDMSG "ERROR: User mask size params are incorrect for source %d\n"

#define dsTGMREGIONSTRINGERR (dsPTINSTRUMENTERROFFSET - 917)
#define dsTGMREGIONSTRINGSEV dsERRSEVFATAL
#define dsTGMREGIONSTRINGSTDMSG "ERROR: Output region string is empty\n"
