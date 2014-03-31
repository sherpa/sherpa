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
#define _DSERROR_PTANALYSIS_H

#define dsPTANALYSISERROFFSET     -4000
#define dsPTANALYSISNUMERRORS     96

/* WRECON ERRORS */
#define dsWRFLUXAFFTERR (dsPTANALYSISERROFFSET - 50)
#define dsWRFLUXAFFTSEV dsERRSEVFATAL
#define dsWRFLUXAFFTSTDMSG "ERROR: FFTConvolve failed in flux_compute() A: %d\n"

#define dsWRFLUXBFFTERR (dsPTANALYSISERROFFSET - 51)
#define dsWRFLUXBFFTSEV dsERRSEVFATAL
#define dsWRFLUXBFFTSTDMSG "ERROR: FFTConvolve failed in flux_compute() B: %d\n"

#define dsWRBADEMAPAXESERR (dsPTANALYSISERROFFSET - 52)
#define dsWRBADEMAPAXESSEV dsERRSEVFATAL
#define dsWRBADEMAPAXESSTDMSG "ERROR: Exposure map has %d != 2 axes\n"

#define dsWRBADEMAPAXLENSERR (dsPTANALYSISERROFFSET - 53)
#define dsWRBADEMAPAXLENSSEV dsERRSEVFATAL
#define dsWRBADEMAPAXLENSSTDMSG "ERROR: Exposure map size (%d, %d) (%s) does not match data size (%d, %d) (%s)\n"

#define dsWRBADFMAPAXESERR (dsPTANALYSISERROFFSET - 54)
#define dsWRBADFMAPAXESSEV dsERRSEVFATAL
#define dsWRBADFMAPAXESSTDMSG "ERROR: Flux image has %d != 2 axes\n"

#define dsWRBADFMAPAXLENSERR (dsPTANALYSISERROFFSET - 55)
#define dsWRBADFMAPAXLENSSEV dsERRSEVFATAL
#define dsWRBADFMAPAXLENSSTDMSG "ERROR: Flux image size (%d, %d) does not match data size (%d, %d)\n"

#define dsWRBADDATAAXESERR (dsPTANALYSISERROFFSET - 56)
#define dsWRBADDATAAXESSEV dsERRSEVFATAL
#define dsWRBADDATAAXESSTDMSG "ERROR: Input data has %d != 2 axes\n"

#define dsWRBADBKGAXESERR (dsPTANALYSISERROFFSET - 57)
#define dsWRBADBKGAXESSEV dsERRSEVFATAL
#define dsWRBADBKGAXESSTDMSG "ERROR: Background has %d != 2 axes\n"

#define dsWRBADBKGAXLENSERR (dsPTANALYSISERROFFSET - 58)
#define dsWRBADBKGAXLENSSEV dsERRSEVFATAL
#define dsWRBADBKGAXLENSSTDMSG "ERROR: Background size (%d, %d) does not match data size (%d, %d)\n"

#define dsWRBADBKGEAXESERR (dsPTANALYSISERROFFSET - 59)
#define dsWRBADBKGEAXESSEV dsERRSEVFATAL
#define dsWRBADBKGEAXESSTDMSG "ERROR: Background error has %d != 2 axes\n"

#define dsWRBADBKGEAXLENSERR (dsPTANALYSISERROFFSET - 60)
#define dsWRBADBKGEAXLENSSEV dsERRSEVFATAL
#define dsWRBADBKGEAXLENSSTDMSG "ERROR: Background error size (%d, %d) does not match data size (%d, %d)\n"

#define dsWRBADCOREAXESERR (dsPTANALYSISERROFFSET - 61)
#define dsWRBADCOREAXESSEV dsERRSEVFATAL
#define dsWRBADCOREAXESSTDMSG "ERROR: Correlation has %d != 2 axes\n"

#define dsWRBADCOREAXLENSERR (dsPTANALYSISERROFFSET - 62)
#define dsWRBADCOREAXLENSSEV dsERRSEVFATAL
#define dsWRBADCOREAXLENSSTDMSG "ERROR: Correlation size (%d, %d) does not match data size (%d, %d)\n"

#define dsWRBADPARAMXSCALESERR (dsPTANALYSISERROFFSET - 63)
#define dsWRBADPARAMXSCALESSEV dsERRSEVFATAL
#define dsWRBADPARAMXSCALESSTDMSG "ERROR: Too many wavelet X scales\n"

#define dsWRBADPARAMYSCALESERR (dsPTANALYSISERROFFSET - 64)
#define dsWRBADPARAMYSCALESSEV dsERRSEVFATAL
#define dsWRBADPARAMYSCALESSTDMSG "ERROR: Too many wavelet Y scales\n"

#define dsWRBADPARAMXYSCALESERR (dsPTANALYSISERROFFSET - 65)
#define dsWRBADPARAMXYSCALESSEV dsERRSEVFATAL
#define dsWRBADPARAMXYSCALESSTDMSG "ERROR: X and y wavelet scale counts do not match\n"

#define dsWROFFSETERR (dsPTANALYSISERROFFSET - 66)
#define dsWROFFSETSEV dsERRSEVFATAL
#define dsWROFFSETSTDMSG "ERROR: Xoffset and yoffset must both be INDEF\n"

#define dsWRCELLTABLEERR (dsPTANALYSISERROFFSET - 67)
#define dsWRCELLTABLESEV dsERRSEVFATAL
#define dsWRCELLTABLESTDMSG "ERROR: Failed to initialize cell table\n"

#define dsWRPSFLIMITERR (dsPTANALYSISERROFFSET - 68)
#define dsWRPSFLIMITSEV dsERRSEVWARNING
#define dsWRPSFLIMITSTDMSG "WARNING: PSF size values not valid beyond %f arcmin off axis\n"

#define dsWRSTKBADNAMEERR (dsPTANALYSISERROFFSET - 69)
#define dsWRSTKBADNAMESEV dsERRSEVFATAL
#define dsWRSTKBADNAMESTDMSG "ERROR: Bad file name in %s stack on line %d\n"

#define dsWRNOGRATERR (dsPTANALYSISERROFFSET - 70)
#define dsWRNOGRATSEV dsERRSEVFATAL
#define dsWRNOGRATSTDMSG "ERROR: Can't read %s keyword in input file.  DM error: %s"

#define dsWRSRCFILECANTERR (dsPTANALYSISERROFFSET - 71)
#define dsWRSRCFILECANTSEV dsERRSEVFATAL
#define dsWRSRCFILECANTSTDMSG "ERROR: Can't create source list file.  DM error: %s"

#define dsWRCOORDCREATEERR (dsPTANALYSISERROFFSET - 72)
#define dsWRCOORDCREATESEV dsERRSEVFATAL
#define dsWRCOORDCREATESTDMSG "ERROR: Could not create coordinate system.  DM error: %s"

#define dsWRCOORDCALCERR (dsPTANALYSISERROFFSET - 73)
#define dsWRCOORDCALCSEV dsERRSEVFATAL
#define dsWRCOORDCALCSTDMSG "ERROR: Could not calculate coordinates.  DM error: %s  Point: (%f, %f)"

#define dsWRCOORDERRCALCERR (dsPTANALYSISERROFFSET - 74)
#define dsWRCOORDERRCALCSEV dsERRSEVFATAL
#define dsWRCOORDERRCALCSTDMSG "ERROR: Could not calculate coordinate errors.  DM error: %s  Point: (%f, %f)"

#define dsWRNOBKGERR (dsPTANALYSISERROFFSET - 75)
#define dsWRNOBKGSEV dsERRSEVFATAL
#define dsWRNOBKGSTDMSG "ERROR: Neither input background nor output background specified"

/* WTRANSFORM ERRORS */
#define dsWTBGFFTERR (dsPTANALYSISERROFFSET - 100)
#define dsWTBGFFTSEV dsERRSEVFATAL
#define dsWTBGFFTSTDMSG "ERROR: FFTConvolve failed in bkgarea() A: %d\n"

#define dsWTCORAFFTERR (dsPTANALYSISERROFFSET - 101)
#define dsWTCORAFFTSEV dsERRSEVFATAL
#define dsWTCORAFFTSTDMSG "ERROR: FFTConvolve failed in correlation() A: %d\n"

#define dsWTCORBFFTERR (dsPTANALYSISERROFFSET - 102)
#define dsWTCORBFFTSEV dsERRSEVFATAL
#define dsWTCORBFFTSTDMSG "ERROR: FFTConvolve failed in correlation() B: %d\n"

#define dsWTCORBGFFTERR (dsPTANALYSISERROFFSET - 103)
#define dsWTCORBGFFTSEV dsERRSEVFATAL
#define dsWTCORBGFFTSTDMSG "ERROR: FFTConvolve failed in background(): %d\n"

#define dsWTCORERRFFTERR (dsPTANALYSISERROFFSET - 104)
#define dsWTCORERRFFTSEV dsERRSEVFATAL
#define dsWTCORERRFFTSTDMSG "ERROR: FFTConvolve failed in correlation_error(): %d\n"

#define dsWTCORERRCORAFFTERR (dsPTANALYSISERROFFSET - 105)
#define dsWTCORERRCORAFFTSEV dsERRSEVFATAL
#define dsWTCORERRCORAFFTSTDMSG "ERROR: FFTConvolve failed in correlation_error_correction() A: %d\n"

#define dsWTCORERRCORBFFTERR (dsPTANALYSISERROFFSET - 106)
#define dsWTCORERRCORBFFTSEV dsERRSEVFATAL
#define dsWTCORERRCORBFFTSTDMSG "ERROR: FFTConvolve failed in correlation_error_correction() B: %d\n"

#define dsWTCORERRBKGERRFFTERR (dsPTANALYSISERROFFSET - 107)
#define dsWTCORERRBKGERRFFTSEV dsERRSEVFATAL
#define dsWTCORERRBKGERRFFTSTDMSG "ERROR: FFTConvolve failed in background_error(): %d\n"

#define dsWTEXPCORAFFTERR (dsPTANALYSISERROFFSET - 108)
#define dsWTEXPCORAFFTSEV dsERRSEVFATAL
#define dsWTEXPCORAFFTSTDMSG "ERROR: FFTConvolve failed in exposure_correction_subr() A: %d\n"

#define dsWTAREACORFFTERR (dsPTANALYSISERROFFSET - 109)
#define dsWTAREACORFFTSEV dsERRSEVFATAL
#define dsWTAREACORFFTSTDMSG "ERROR: FFTConvolve failed in area_correction(): %d\n"

#define dsWTBADPARAMITERERR (dsPTANALYSISERROFFSET - 110)
#define dsWTBADPARAMITERSEV dsERRSEVFATAL
#define dsWTBADPARAMITERSTDMSG "ERROR: Bad maxiter parameter value (%d)\n"

#define dsWTBADPARAMSTOPERR (dsPTANALYSISERROFFSET - 111)
#define dsWTBADPARAMSTOPSEV dsERRSEVFATAL
#define dsWTBADPARAMSTOPSTDMSG "ERROR: Bad iterstop parameter value (%f)\n"

#define dsWTBADPARAMSIGERR (dsPTANALYSISERROFFSET - 112)
#define dsWTBADPARAMSIGSEV dsERRSEVFATAL
#define dsWTBADPARAMSIGSTDMSG "ERROR: Bad sigthresh parameter value (%f)\n"

#define dsWTBADPARAMSIGBNDSERR (dsPTANALYSISERROFFSET - 113)
#define dsWTBADPARAMSIGBNDSSEV dsERRSEVFATAL
#define dsWTBADPARAMSIGBNDSSTDMSG "ERROR: Source significance is outside supported bounds - spurious results may occur.\n"

#define dsWTBADPARAMBKGSIGERR (dsPTANALYSISERROFFSET - 114)
#define dsWTBADPARAMBKGSIGSEV dsERRSEVFATAL
#define dsWTBADPARAMBKGSIGSTDMSG "ERROR: Bad bkgsigthresh parameter value (%f)\n"

#define dsWTBADPARAMBKGSIGBNDSERR (dsPTANALYSISERROFFSET - 115)
#define dsWTBADPARAMBKGSIGBNDSSEV dsERRSEVFATAL
#define dsWTBADPARAMBKGSIGBNDSSTDMSG "ERROR: Background significance is outside supported bounds - spurious results may occur.\n"

#define dsWTCREATECMERR (dsPTANALYSISERROFFSET - 116)
#define dsWTCREATECMSEV dsERRSEVFATAL
#define dsWTCREATECMSTDMSG "ERROR: Couldn't create correlation maxima file.  DM error: %s"

#define dsWTEXPCORBFFTERR (dsPTANALYSISERROFFSET - 117)
#define dsWTEXPCORBFFTSEV dsERRSEVFATAL
#define dsWTEXPCORBFFTSTDMSG "ERROR: FFTConvolve failed in exposure_correction_subr() B: %d\n"

#define dsWTBADPARAMXSCALESERR (dsPTANALYSISERROFFSET - 118)
#define dsWTBADPARAMXSCALESSEV dsERRSEVFATAL
#define dsWTBADPARAMXSCALESSTDMSG "ERROR: Too many X wavelet scales\n"

#define dsWTBADPARAMYSCALESERR (dsPTANALYSISERROFFSET - 119)
#define dsWTBADPARAMYSCALESSEV dsERRSEVFATAL
#define dsWTBADPARAMYSCALESSTDMSG "ERROR: Too many Y wavelet scales\n"

#define dsWTBADPARAMXYSCALESERR (dsPTANALYSISERROFFSET - 120)
#define dsWTBADPARAMXYSCALESSEV dsERRSEVFATAL
#define dsWTBADPARAMXYSCALESSTDMSG "ERROR: X and y wavelet scale counts do not match\n"


#define dsWTBADEMAPAXESERR (dsPTANALYSISERROFFSET - 121)
#define dsWTBADEMAPAXESSEV dsERRSEVFATAL
#define dsWTBADEMAPAXESSTDMSG "ERROR: exposure map has %d != 2 axes\n"

#define dsWTBADEMAPAXLENSERR (dsPTANALYSISERROFFSET - 122)
#define dsWTBADEMAPAXLENSSEV dsERRSEVFATAL
#define dsWTBADEMAPAXLENSSTDMSG "ERROR: Exposure map size (%d, %d) (%s) does not match data size (%d, %d) (%s)\n"

#define dsWTBADPARAMWSCALESIZEERR (dsPTANALYSISERROFFSET - 123)
#define dsWTBADPARAMWSCALESIZESEV dsERRSEVFATAL
#define dsWTBADPARAMWSCALESIZESTDMSG "ERROR: maximum wavelet scale (%f pixels) too big for dataset size (%d, %d)\n"

/* CELLDETECT ERRORS */
#define dsCDBADBGAXESERR (dsPTANALYSISERROFFSET - 150)
#define dsCDBADBGAXESSEV dsERRSEVFATAL
#define dsCDBADBGAXESSTDMSG "ERROR: Background file has %d != 2 axes\n"

#define dsCDBADBGAXLENSERR (dsPTANALYSISERROFFSET - 151)
#define dsCDBADBGAXLENSSEV dsERRSEVFATAL
#define dsCDBADBGAXLENSSTDMSG "ERROR: Background size (%d, %d) does not match data size (%d, %d)\n"

#define dsCDWCSINSTRERR (dsPTANALYSISERROFFSET - 152)
#define dsCDWCSINSTRSEV dsERRSEVFATAL

#define dsCDWCSINSTRSTDMSG "ERROR: WCS information or TELESCOP/INSTRUMEN/DETNAM not present in input: cannot run variable kernel\n"

#define dsCDPSFCANTERR (dsPTANALYSISERROFFSET - 153)
#define dsCDPSFCANTSEV dsERRSEVWARNING
#define dsCDPSFCANTSTDMSG "WARNING: Failed to calculate psf size for off axis angle %f\n"

#define dsCDCELLSIZEERR (dsPTANALYSISERROFFSET - 154)
#define dsCDCELLSIZESEV dsERRSEVFATAL
#define dsCDCELLSIZESTDMSG "ERROR: Fixed cell size must be 1 or a multiple or 3!\n"

#define dsCDOFFSETERR (dsPTANALYSISERROFFSET - 155)
#define dsCDOFFSETSEV dsERRSEVFATAL
#define dsCDOFFSETSTDMSG "ERROR: Xoffset and yoffset must both be INDEF\n"

#define dsCDCELLTABLEERR (dsPTANALYSISERROFFSET - 156)
#define dsCDCELLTABLESEV dsERRSEVFATAL
#define dsCDCELLTABLESTDMSG "ERROR: Failed to initialize cell table\n"

#define dsCDPSFLIMITERR (dsPTANALYSISERROFFSET - 157)
#define dsCDPSFLIMITSEV dsERRSEVWARNING
#define dsCDPSFLIMITSTDMSG "WARNING: PSF size values not valid beyond %f arcmin off axis\n"

#define dsCDBADINAXESERR (dsPTANALYSISERROFFSET - 158)
#define dsCDBADINAXESSEV dsERRSEVFATAL
#define dsCDBADINAXESSTDMSG "ERROR: Input file has %d != 2 axes\n"

#define dsCDCOORDCREATEERR (dsPTANALYSISERROFFSET - 159)
#define dsCDCOORDCREATESEV dsERRSEVFATAL
#define dsCDCOORDCREATESTDMSG "ERROR: Could not create coordinate system.  DM error: %s"

#define dsCDCOORDCALCERR (dsPTANALYSISERROFFSET - 160)
#define dsCDCOORDCALCSEV dsERRSEVFATAL
#define dsCDCOORDCALCSTDMSG "ERROR: Could not calculate coordinates.  DM error: %s  Point: (%f, %f)"

#define dsCDCOORDERRCALCERR (dsPTANALYSISERROFFSET - 161)
#define dsCDCOORDERRCALCSEV dsERRSEVFATAL
#define dsCDCOORDERRCALCSTDMSG "ERROR: Could not calculate coordinate errors.  DM error: %s  Point: (%f, %f)"

#define dsCDNOEXPMAPERR (dsPTANALYSISERROFFSET - 162)
#define dsCDNOEXPMAPSEV dsERRSEVFATAL
#define dsCDNOEXPMAPSTDMSG "ERROR: no exposure map specifed"

#define dsCDFEWEXPMAPERR (dsPTANALYSISERROFFSET - 163)
#define dsCDFEWEXPMAPSEV dsERRSEVFATAL
#define dsCDFEWEXPMAPSTDMSG "ERROR: too few exposure maps in list"

#define dsCDSIZEEXPMAPERR (dsPTANALYSISERROFFSET - 164)
#define dsCDSIZEEXPMAPSEV dsERRSEVFATAL
#define dsCDSIZEEXPMAPSTDMSG "ERROR: exposure map size (%ld,%ld) (%s) does not match image size (%ld,%ld) (%s)"



/* STK_READ_NUM ERRORS */
#define dsSTKREADNUMLOWNUMERR (dsPTANALYSISERROFFSET - 200 )
#define dsSTKREADNUMLOWNUMSEV dsERRSEVWARNING
#define dsSTKREADNUMLOWNUMSTDMSG "WARNING: Stack element < 1, reseting num to 1.\n"

#define dsSTKREADNUMHIGHNUMERR (dsPTANALYSISERROFFSET - 201 )
#define dsSTKREADNUMHIGHNUMSEV dsERRSEVFATAL
#define dsSTKREADNUMHIGHNUMSTDMSG "ERROR: Number specified is greater than stack count.\n"

#define dsSTKREADNUMEMPTYERR (dsPTANALYSISERROFFSET - 202 )
#define dsSTKREADNUMEMPTYSEV dsERRSEVWARNING
#define dsSTKREADNUMEMPTYSTDMSG "WARNING: %ith element of stack is empty.\n"



/* BKGD CALC GLOBAL COUNT RATE SPECIFIC ERRORS */
#define dsFITPARAMERR (dsPTANALYSISERROFFSET - 400)
#define dsFITPARAMSEV dsERRSEVFATAL
#define dsFITPARAMSTDMSG "ERROR: Fitting parameter error.\n"

#define dsHISTLISTERR (dsPTANALYSISERROFFSET - 401)
#define dsHISTLISTSEV dsERRSEVFATAL
#define dsHISTLISTSTDMSG "ERROR: Failed to create the histogram linked list.\n"

#define dsWRITEBACKGDERR (dsPTANALYSISERROFFSET - 402)
#define dsWRITEBACKGDSEV dsERRSEVFATAL
#define dsWRITEBACKGDSTDMSG "ERROR: Failed to write the background events.\n"


/* TCD SPECIFIC ERRORS */
#define dsTCDNAXESERR (dsPTANALYSISERROFFSET - 650)
#define dsTCDNAXESSEV dsERRSEVFATAL
#define dsTCDNAXESSTDMSG "ERROR: User passed nmuber of axes <= 0.\n"

#define dsTCDPADLTOLDERR (dsPTANALYSISERROFFSET -651)
#define dsTCDPADLTOLDSEV dsERRSEVFATAL
#define dsTCDPADLTOLDSTDMSG "ERROR: Required padding less than original image size.\n"

#define dsTCDUNKWNPADERR (dsPTANALYSISERROFFSET - 652)
#define dsTCDUNKWNPADSEV dsERRSEVFATAL
#define dsTCDUNKWNPADSTDMSG "ERROR: Unknown padding specified.\n"

#define dsTCDLAXES0ERR (dsPTANALYSISERROFFSET - 653)
#define dsTCDLAXES0SEV dsERRSEVFATAL
#define dsTCDLAXES0STDMSG "ERROR: User specified lAxes[i] = 0.\n"

#define dsTCDUNSUPPORTNAXESERR (dsPTANALYSISERROFFSET - 654)
#define dsTCDUNSUPPORTNAXESSEV dsERRSEVFATAL
#define dsTCDUNSUPPORTNAXESSTDMSG "ERROR: Kernel lib doesn't support specified nAxes.\n"

#define dsTCDNAXESMISMATCHERR (dsPTANALYSISERROFFSET - 655)
#define dsTCDNAXESMISMATCHSEV dsERRSEVFATAL
#define dsTCDNAXESMISMATCHSTDMSG "ERROR: Kernel and image nAxes don't match.\n"

#define dsTCDINCONSISTENTERR (dsPTANALYSISERROFFSET - 656)
#define dsTCDINCONSISTENTSEV dsERRSEVFATAL
#define dsTCDINCONSISTENTSTDMSG "ERROR: Length of axes specified in string not same.\n"

#define dsTCDUNKWNKERNELERR (dsPTANALYSISERROFFSET - 657)
#define dsTCDUNKWNKERNELSEV dsERRSEVFATAL
#define dsTCDUNKWNKERNELSTDMSG "ERROR: Unknown TCD kernel type: %s.\n"


/* CSMOOTH ERRORS */
#define dsSIZEMISMATCHERR (dsPTANALYSISERROFFSET - 850)
#define dsSIZEMISMATCHSEV dsERRSEVFATAL
#define dsSIZEMISMATCHSTDMSG "ERROR: The size of input files %s is different from the size of image file %s.\n"

#define dsBKGMODEERR (dsPTANALYSISERROFFSET - 851)
#define dsBKGMODESEV dsERRSEVWARNING
#define dsBKGMODESTDMSG "WARNING: Bkgmode is set to user, bgrdcode cannot be -1.\n"


/* MKRMF ERRORS */
#define dsMKRMFDIMENSIONERR (dsPTANALYSISERROFFSET - 900)
#define dsMKRMFDIMENSIONSEV  dsERRSEVFATAL
#define dsMKRMFDIMENSIONMSG "ERROR: Dimension (%d) exceeds the current maximum capability (%d).\n"

#define dsMKRMFINFILESYNTXERR (dsPTANALYSISERROFFSET - 901)
#define dsMKRMFINFILESYNTXSEV  dsERRSEVFATAL
#define dsMKRMFINFILESYNTXMSG "ERROR: Invalid syntax for input file, \"%s\".\n"

#define dsMKRMFUSRGRIDSYNTXERR (dsPTANALYSISERROFFSET - 902)
#define dsMKRMFUSRGRIDSYNTXSEV  dsERRSEVFATAL
#define dsMKRMFUSRGRIDSYNTXMSG "ERROR: Invalid syntax for grids input, \"%s\".\n\t Suggest syntax: %s.\n"

#define dsMKRMFMISMTCHNMESERR (dsPTANALYSISERROFFSET - 903)
#define dsMKRMFMISMTCHNMESSEV  dsERRSEVFATAL
#define dsMKRMFMISMTCHNMESMSG "ERROR: Either or all of axis names (%s) do not match column names (%s) in \"%s\" file.\n"

#define dsMKRMFOUTPUTFMTERR (dsPTANALYSISERROFFSET - 904)
#define dsMKRMFOUTPUTFMTSEV  dsERRSEVFATAL
#define dsMKRMFOUTPUTFMTMSG "ERROR: Invalid RMF output format, \"%s\".\n Suggested value: \"legacy\" or \"cxc\".\n"

#define dsMKRMFNOTLOGBINERR (dsPTANALYSISERROFFSET - 905)
#define dsMKRMFNOTLOGBINSEV  dsERRSEVFATAL
#define dsMKRMFNOTLOGBINMSG "ERROR: Not supported logarithmic bin for \"%s\" axis.\n"

#define dsMKRMFZEROBINBUMERR (dsPTANALYSISERROFFSET - 906)
#define dsMKRMFZEROBINBUMSEV  dsERRSEVFATAL
#define dsMKRMFZEROBINBUMMSG "ERROR: Zero bin number for \"%s\".\n\t DM Error: %s.\n"

#define dsMKRMFDECBINSTPERR (dsPTANALYSISERROFFSET - 907)
#define dsMKRMFDECBINSTPSEV  dsERRSEVFATAL
#define dsMKRMFDECBINSTPMSG "ERROR: Do not allow decimal number for bin step (%f) in \"%s\".\n"

#define dsMKRMFNEGVALUEERR (dsPTANALYSISERROFFSET - 908)
#define dsMKRMFNEGVALUESEV  dsERRSEVFATAL
#define dsMKRMFNEGVALUEMSG "ERROR: Negative value calculated for \"%s\".\n"

#define dsMKRMFLONGKWDERR (dsPTANALYSISERROFFSET - 909)
#define dsMKRMFLONGKWDSEV  dsERRSEVWARNING
#define dsMKRMFLONGKWDMSG "WARNING: Keyword, \"%s\", should limite to %d characters.\n"

#define dsMKRMFGENERICERR (dsPTANALYSISERROFFSET - 910)
#define dsMKRMFGENERICSEV  dsERRSEVWARNING
#define dsMKRMFGENERICMSG "WARNING: Unspecified mkrmf error.\n"

#define dsMKRMFINFILEREGNUMERR (dsPTANALYSISERROFFSET - 911)
#define dsMKRMFINFILEREGNUMSEV  dsERRSEVFATAL
#define dsMKRMFINFILEREGNUMMSG "ERROR: Ambiguous data filter for input file, \"%s\"\n Suggestion:\n\tspecify properate data filter(s) if none, or\n\tvalidate, for example,  REGNUM value if any.\n"


/* LIGHTCURVE ERRORS */
#define dsBKGGDBINERR (dsPTANALYSISERROFFSET - 950)
#define dsBKGGDBINSEV dsERRSEVFATAL
#define dsBKGGDBINSTDMSG "ERROR: Problem binning background counts.\n"

#define dsSETTIMEBINSERR (dsPTANALYSISERROFFSET - 951)
#define dsSETTIMEBINSSEV dsERRSEVFATAL
#define dsSETTIMEBINSSTDMSG "ERROR: Failed to set timing bins.\n"

#define dsSOURCEBINERR (dsPTANALYSISERROFFSET - 952)
#define dsSOURCEBINSEV dsERRSEVFATAL
#define dsSOURCEBINSTDMSG "ERROR: Failed to bin source counts.\n"


