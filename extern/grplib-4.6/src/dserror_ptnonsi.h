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
#define _DSERROR_PTNONSI_H

#define dsPTNONSIERROFFSET     -5000
#define dsPTNONSINUMERRORS     61

/* MTL_BUILD_TABLE ERRORS */

#define dsMTLLOOKUPCOLERR (dsPTNONSIERROFFSET - 0)
#define dsMTLLOOKUPCOLSEV dsERRSEVFATAL
#define dsMTLLOOKUPCOLSTDMSG "ERROR: Lookup table does not contain required columns.\n"

#define dsMTLSEEDFILEERR (dsPTNONSIERROFFSET - 1)
#define dsMTLSEEDFILESEV dsERRSEVFATAL
#define dsMTLSEEDFILESTDMSG "ERROR: No seed files specified.\n"

#define dsMTLDATALOOKUPERR (dsPTNONSIERROFFSET - 2)
#define dsMTLDATALOOKUPSEV dsERRSEVWARNING
#define dsMTLDATALOOKUPSTDMSG "WARNING: The lookup table specified contains no data for the file %s.\n"

#define dsMTLDATALIMERR (dsPTNONSIERROFFSET - 3)
#define dsMTLDATALIMSEV dsERRSEVWARNING
#define dsMTLDATALIMSTDMSG "WARNING: Specified obistart and stop are not contained within the data provided.  Using %f to %f.\n"


/* MTL_BUILD_GTI ERRORS */
#define dsMBGSMOOTHERR (dsPTNONSIERROFFSET - 50)
#define dsMBGSMOOTHSEV dsERRSEVWARNING
#define dsMBGSMOOTHSTDMSG "WARNING: Smooth factor for column %s was set higher than one, and it is not an integer or real column.  Setting smooth factor to one.\n"

#define dsMBGLIMITCONFLICTERR (dsPTNONSIERROFFSET - 51)
#define dsMBGLIMITCONFLICTSEV dsERRSEVWARNING
#define dsMBGLIMITCONFLICTSTDMSG "WARNING: there has been some sort of data type conflict in the limits.\nThe remaining limits will be re-evaluated.  This may result in additional limits conditions being discarded.\nThe following errors occurred:\n"

#define dsMBGLIMITSERR (dsPTNONSIERROFFSET - 52)
#define dsMBGLIMITSSEV dsERRSEVWARNING
#define dsMBGLIMITSSTDMSG "WARNING: the following errors occurred during initialization of the limits.  They are not fatal errors - the tool will continue.\n"


/* ASP_CALC_OFFSETS */
#define dsACOGAPWERR (dsPTNONSIERROFFSET - 100)
#define dsACOGAPWSEV dsERRSEVWARNING
#define dsACOGAPWSTDMSG "WARNING: Gap in the data from time %f to time %f.\n"

#define dsACODETNAMEERR (dsPTNONSIERROFFSET - 101)
#define dsACODETNAMESEV dsERRSEVFATAL
#define dsACODETNAMESTDMSG "ERROR: The DETNAM %s, does not match HRC-I, HRC-S or any ACIS string. Please check that DETNAM is set correctly.\n"

#define dsACONODELERR (dsPTNONSIERROFFSET - 102)
#define dsACONODELSEV dsERRSEVFATAL
#define dsACONODELSTDMSG "ERROR: Failed to find the time between records in the aspect solution file.\n"

#define dsACONOPARAMERR (dsPTNONSIERROFFSET - 103)
#define dsACONOPARAMSEV dsERRSEVFATAL
#define dsACONOPARAMSTDMSG "ERROR: One or more parameters not found in %s.\n"


/* ASPHIST ERRORS */
#define dsASPHISTNEGERR (dsPTNONSIERROFFSET - 150)
#define dsASPHISTNEGSEV dsERRSEVWARNING
#define dsASPHISTNEGSTDMSG "WARNING: Negative value illegal for octree: %f\n"

#define dsASPHISTLTFERR (dsPTNONSIERROFFSET - 151)
#define dsASPHISTLTFSEV dsERRSEVFATAL
#define dsASPHISTLTFSTDMSG "ERROR: Problem with livetime correction file %s\n"


/* EPHEM_CALC_VG_ANGLES ERRORS*/
#define dsECVACOSERR (dsPTNONSIERROFFSET - 200)
#define dsECVACOSSEV dsERRSEVFATAL
#define dsECVACOSSTDMSG "ERROR: Tried to take arccosine of a number not between -1 and 1.\n"

#define dsECVARUNOUTWERR (dsPTNONSIERROFFSET - 201)
#define dsECVARUNOUTWSEV dsERRSEVWARNING
#define dsECVARUNOUTWSTDMSG "WARNING: Files ended before requested stop time. Output file %s was generated, but ends at time %f.\n"


/* TEL_DATA_CLEAN ERRORS */
#define dsTDCINPUTSTKERR (dsPTNONSIERROFFSET - 300)
#define dsTDCINPUTSTKSEV dsERRSEVFATAL
#define dsTDCINPUTSTKSTDMSG "ERROR: Number of input data files has to be equal to the number of input status files. \n"


/* SIM_DATA_SELECT ERRORS */
#define dsSIM2UPDATEDERR (dsPTNONSIERROFFSET - 350)
#define dsSIM2UPDATEDSEV dsERRSEVWARNING
#define dsSIM2UPDATEDSTDMSG "ERROR: Both samples are UPDATED. Please check input file."

#define dsSIMNOUPDATEDERR (dsPTNONSIERROFFSET - 351)
#define dsSIMNOUPDATEDSEV dsERRSEVWARNING
#define dsSIMNOUPDATEDSTDMSG "WARNING: Neither sample was UPDATED. Please check input file."

#define dsSDSNOROWSERR (dsPTNONSIERROFFSET - 352) 
#define dsSDSNOROWSSEV dsERRSEVFATAL
#define dsSDSNOROWSSTDMSG "ERROR: No rows written to output.\n"

#define dsSDSNOINERR (dsPTNONSIERROFFSET - 353) 
#define dsSDSNOINSEV dsERRSEVFATAL
#define dsSDSNOINSTDMSG "ERROR: No input files to process.\n"

#define dsSDSNEEDROWERR (dsPTNONSIERROFFSET - 354)
#define dsSDSNEEDROWSEV  dsERRSEVWARNING
#define dsSDSNEEDROWSTDMSG "WARNING: Not enough rows to process input file list:\n %s.\n Possible Causes: Bad Quality flags in the input files\n\t\t  Time filter does not include times in files.\n\t\t  3 or Less rows in input files.\n"

#define dsSDSBADMATCHERR (dsPTNONSIERROFFSET - 355)
#define dsSDSBADMATCHSEV dsERRSEVFATAL
#define dsSDSBADMATCHSTDMSG "ERROR: The time and major frame numbers differences do not correlate in the input file at time %f vs %f"

#define dsSDSPREVMJFERR (dsPTNONSIERROFFSET - 356)
#define dsSDSPREVMJFSEV dsERRSEVWARNING
#define dsSDSPREVMJFSTDMSG "Warning: Major frame does not match the previous. Cannot interpolate. Not Updating.\n"

#define dsSDSPREVTIMEERR (dsPTNONSIERROFFSET - 357)
#define dsSDSPREVTIMESEV dsERRSEVWARNING
#define dsSDSPREVTIMESTDMSG "Warning: Time does not match the previous. Cannot interpolate. Not Updating.\n"


/* SIM_COMPUTE_STF_POS ERRORS */
#define dsSCPTIMESERR (dsPTNONSIERROFFSET - 500)
#define dsSCPTIMESSEV dsERRSEVFATAL
#define dsSCPTIMESSTDMSG "ERROR: The intersected times to process on are zeroes.\n\t Please check input data.\n"

#define dsSCPSTATUSERR (dsPTNONSIERROFFSET - 501)
#define dsSCPSTATUSSEV dsERRSEVFATAL
#define dsSCPSTATUSSTDMSG "ERROR: NULL pointers passed to check_status function.\n"

#define dsSCPBADCOLERR (dsPTNONSIERROFFSET - 502)
#define dsSCPBADCOLSEV dsERRSEVFATAL
#define dsSCPBADCOLSTDMSG "ERROR: Bad column Number passed back.\n"


/* OBI_MERGE_TIMES ERRORS */
#define dsOMTBLOBJCOLERR (dsPTNONSIERROFFSET - 550)
#define dsOMTBLOBJCOLSEV dsERRSEVWARNING
#define dsOMTBLOBJCOLSTDMSG "WARNING: %s is not a column in the specified obiBlock object.\n"

#define dsOMTNULLSSERR (dsPTNONSIERROFFSET - 551)
#define dsOMTNULLSSSEV dsERRSEVFATAL
#define dsOMTNULLSSSTDMSG "ERROR: No time subspaces found, null intersection of inputs.\n"


/* OBI_SIM_LOOK ERRORS */
#define dsOSLBADRNGERR (dsPTNONSIERROFFSET - 600)
#define dsOSLBADRNGSEV dsERRSEVWARNING
#define dsOSLBADRNGSTDMSG "WARNING: SIM position is not within expected boundaries\n"

/* OBI_ENG_LOOK ERRORS */
#define dsOELNOIDERR (dsPTNONSIERROFFSET - 650)
#define dsOELNOIDSEV dsERRSEVFATAL
#define dsOELNOIDSTDMSG "ERROR: Couldn't find obsid %d bewteen %f and %f.\n"

#define dsOELNOALTERR (dsPTNONSIERROFFSET - 651) 
#define dsOELNOALTSEV dsERRSEVWARNING
#define dsOELNOALTSTDMSG "WARNING: Angle files stack has emptied at time %f.\n"

#define dsOELFAULTWERR (dsPTNONSIERROFFSET - 652)
#define dsOELFAULTWSEV dsERRSEVWARNING
#define dsOELFAULTWSTDMSG "WARNING: Grating value is FAULT.\n"

#define dsOELUKNWERR (dsPTNONSIERROFFSET - 653)
#define dsOELUKNWSEV dsERRSEVWARNING
#define dsOELUKNWSTDMSG "WARNING: Possible data gap. Flagged as bad status.\n"


/* OBI_ASP_LOOK ERRORS */
#define dsOASPSKIPWERR (dsPTNONSIERROFFSET - 700)
#define dsOASPSKIPWSEV dsERRSEVWARNING
#define dsOASPSKIPWSTDMSG "WARNING: Skipping input row with time length shorter than %f.\n"

#define dsOASPNOIDERR (dsPTNONSIERROFFSET - 701)
#define dsOASPNOIDSEV dsERRSEVFATAL
#define dsOASPNOIDSTDMSG "ERROR: Requested obsid %d was not found in the input files.\n"

#define dsOASPNOKALWERR (dsPTNONSIERROFFSET - 702)
#define dsOASPNOKALWSEV dsERRSEVWARNING
#define dsOASPNOKALWSTDMSG "WARNING: Requested obsid %d was not found in aspect mode KALMAN in this file.\n"

#define dsOASPNOROWSERR (dsPTNONSIERROFFSET - 703)
#define dsOASPNOROWSSEV dsERRSEVFATAL
#define dsOASPNOROWSSTDMSG "ERROR: Appear to have no rows in input file.\n"

#define dsOASPEARLYWERR (dsPTNONSIERROFFSET - 704)
#define dsOASPEARLYWSEV dsERRSEVWARNING
#define dsOASPEARLYWSTDMSG "WARNING: Passed beyond the tstart. Collecting previous row time.\n\t May need to increase tstart time.\n"

#define dsOASPLATEWERR (dsPTNONSIERROFFSET - 705)
#define dsOASPLATEWSEV dsERRSEVWARNING
#define dsOASPLATEWSTDMSG "WARNING: Passed beyond the tstop. Collecting previous row time.\n\t May need to increase tstop time.\n"


/* OBI_CHECK_TOL ERRORS */
#define dsOCTTOLVALERR (dsPTNONSIERROFFSET - 800)
#define dsOCTTOLVALSEV dsERRSEVWARNING
#define dsOCTTOLVALSTDMSG  "WARNING: Could not find a matching tolerance value for element %s, value %s.  Check tolerance file against input data file.\n"

#define dsOCTTOLTYPEERR (dsPTNONSIERROFFSET - 801)
#define dsOCTTOLTYPESEV dsERRSEVFATAL
#define dsOCTTOLTYPESTDMSG "ERROR: The data type of element %s does not match the input file column type.\n"

#define dsOCTDURATIONERR (dsPTNONSIERROFFSET - 802)
#define dsOCTDURATIONSEV dsERRSEVWARNING
#define dsOCTDURATIONSTDMSG "WARNING:  The obi element from %20.10g to %20.10g does not match\nthe nominal obi, and has a duration greater than %20.10g.\n"

#define dsOCTOUTPARERR (dsPTNONSIERROFFSET - 803)
#define dsOCTOUTPARSEV dsERRSEVFATAL
#define dsOCTOUTPARSTDMSG "ERROR: Input file %s did not contain the necessary parameters, output file %s not created.\n"

#define dsOCTNOMMATCHERR (dsPTNONSIERROFFSET - 804)
#define dsOCTNOMMATCHSEV dsERRSEVWARNING
#define dsOCTNOMMATCHSTDMSG "WARNING: Could not find any obi elements that matched the nominal values.\n"

#define dsOCTBADTOLERR (dsPTNONSIERROFFSET - 805)
#define dsOCTBADTOLSEV dsERRSEVWARNING
#define dsOCTBADTOLSTDMSG "WARNING: Partial or alternate matches are nonsensical for %s data.  Review tolerances file.\n"

#define dsOCTTOLSUPPORTERR (dsPTNONSIERROFFSET - 806)
#define dsOCTTOLSUPPORTSEV dsERRSEVWARNING
#define dsOCTTOLSUPPORTSTDMSG "WARNING: The only check data types that are currently supported are shorts, longs, and character strings.\n\tCheck the tolerances file.\n"

#define dsOCTFINDTOLERR (dsPTNONSIERROFFSET - 807)
#define dsOCTFINDTOLSEV dsERRSEVFATAL
#define dsOCTFINDTOLSTDMSG "ERROR: Could not find element %s in tolerance file.\n"

#define dsOCTCHECKERR (dsPTNONSIERROFFSET - 808)
#define dsOCTCHECKSEV dsERRSEVFATAL
#define dsOCTCHECKSTDMSG   "ERROR: Could not perform tolerance check for element %s.\nTolerance file should have both an exact and a wrong case.\n"


/* CALC_MEAN_POINTING ERRORS */
#define dsCMPNOSTARTERR (dsPTNONSIERROFFSET - 850)
#define dsCMPNOSTARTSEV dsERRSEVFATAL
#define dsCMPNOSTARTSTDMSG "ERROR: Times in AIPROPS file are all before tstart.\n"

#define dsCMPRUNOUTERR (dsPTNONSIERROFFSET - 851)
#define dsCMPRUNOUTSEV dsERRSEVWARNING
#define dsCMPRUNOUTSTDMSG "WARNING: AIPROPS files end at %f which is before\n tstop(%f).\n"

#define dsCMPINKALERR (dsPTNONSIERROFFSET - 852)
#define dsCMPINKALSEV dsERRSEVWARNING
#define dsCMPINKALSTDMSG "WARNING: AIPROPS files end in a KALMAN interval.\n"


/* OBI_ACISHK_LOOK ERRORS */
#define dsAHLINSTERR (dsPTNONSIERROFFSET - 900)
#define dsAHLINSTSEV dsERRSEVFATAL
#define dsAHLINSTSTDMSG "ERROR: Instrument must be ACIS.\n"

#define dsAHLRDMODEERR (dsPTNONSIERROFFSET - 901)
#define dsAHLRDMODESEV dsERRSEVFATAL
#define dsAHLRDMODESTDMSG "ERROR: Cannot identify this readmode: %s.\n"

#define dsAHLNOTIMESWERR (dsPTNONSIERROFFSET - 902)
#define dsAHLNOTIMESWSEV dsERRSEVWARNING
#define dsAHLNOTIMESWSTDMSG "WARNING: File %s does not intersect with given start and stop times.\n"

/* OBI_HRCHK_LOOK ERRORS */
#define dsOHLINSTERR (dsPTNONSIERROFFSET - 950)
#define dsOHLINSTSEV dsERRSEVFATAL
#define dsOHLINSTSTDMSG "ERROR: Instrument must be HRC.\n"

#define dsOHLREADCOLHERR (dsPTNONSIERROFFSET - 951)
#define dsOHLREADCOLHSEV dsERRSEVFATAL
#define dsOHLREADCOLHSTDMSG "ERROR: Couldn't read required columns from the housekeeping file,%s.\n"

#define dsOHLREADCOLEERR (dsPTNONSIERROFFSET - 952)
#define dsOHLREADCOLESEV dsERRSEVFATAL
#define dsOHLREADCOLESTDMSG "ERROR: Couldn't read required columns from the engineering file,%s.\n"

#define dsOHLDETNAMEERR (dsPTNONSIERROFFSET - 953)
#define dsOHLDETNAMESEV dsERRSEVWARNING
#define dsOHLDETNAMESTDMSG "WARNING: Invalid detector %s, resetting to NONE.\n"

#define dsOHLNOTIMEERR (dsPTNONSIERROFFSET - 954)
#define dsOHLNOTIMESEV dsERRSEVWARNING
#define dsOHLNOTIMESTDMSG "WARNING: Did not find any data in the time interval.\n"

/* OBSPAR_UPD ERRORS */
#define dsUPDNODETERR (dsPTNONSIERROFFSET - 1000)
#define dsUPDNODETSEV dsERRSEVWARNING
#define dsUPDNODETSTDMSG "WARNING: detector name is either NONE or not found. Setting sim values to 0.0\n"


