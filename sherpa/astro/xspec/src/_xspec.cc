// 
//  Copyright (C) 2007  Smithsonian Astrophysical Observatory
//
//
//  This program is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License along
//  with this program; if not, write to the Free Software Foundation, Inc.,
//  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
//

// INIT_XSPEC is used by xspecmodelfct() in xspec_extension.hh, so it needs to
// be defined before that file is included
int _sherpa_init_xspec_library();
#define INIT_XSPEC _sherpa_init_xspec_library

#include <iostream>
#include <fstream>
#include "sherpa/astro/xspec_extension.hh"

#define ABUND_SIZE (30) // number of elements in Solar Abundance table

extern "C" {
void init_xspec();

// Lifted from XSPEC 12 include directory
char* FGXSCT(void); 
char* FGSOLR(void);

char* FGMSTR(char* dname);
float FGABND(char* element);

void FPSOLR(const char* table, int* ierr);
void FPXSCT(const char* csection, int* ierr);
void FPMSTR(const char* value1, const char* value2);
void FPSLFL(float rvalue[], int nvalue, int *ierr);

void FNINIT(void);
float csmgq0(void);
float csmgh0(void);
float csmgl0(void);
void csmpq0(float q0);
void csmph0(float H0);
void csmpl0(float lambda0);
int FGCHAT();
void FPCHAT(int chat);
int xs_getVersion(char* buffer, int buffSize);

void xsaped_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xsbape_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xsblbd_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xsbbrd_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
//void xsbexrav_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void C_xsbexrav(const double* energy, int nFlux, const double* params, int spectrumNumber, double* flux, double* fluxError, const char* initStr);
void C_xsbexriv(const double* energy, int nFlux, const double* params, int spectrumNumber, double* flux, double* fluxError, const char* initStr);
void C_brokenPowerLaw(const double* energy, int nFlux, const double* params, int spectrumNumber, double* flux, double* fluxError, const char* initStr);
void C_broken2PowerLaw(const double* energy, int nFlux, const double* params, int spectrumNumber, double* flux, double* fluxError, const char* initStr);
void C_sirf(const double* energy, int nFlux, const double* params, int spectrumNumber, double* flux, double* fluxError, const char* initStr);
void xsbmc_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xsbrms_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xsbvpe_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void c6mekl_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void c6pmekl_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void c6pvmkl_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void c6vmekl_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void cemekl_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void C_cemVMekal(const double* energy, int nFlux, const double* params, int spectrumNumber, double* flux, double* fluxError, const char* initStr);
// void xscflw_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void C_xscflw(const double* energy, int nFlux, const double* params, int spectrumNumber, double* flux, double* fluxError, const char* initStr);
void compbb_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void compls_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
// void xscompps_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void C_xscompps(const double* energy, int nFlux, const double* params, int spectrumNumber, double* flux, double* fluxError, const char* initStr);
void compst_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xstitg_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void C_cutoffPowerLaw(const double* energy, int nFlux, const double* params, int spectrumNumber, double* flux, double* fluxError, const char* initStr);
void disk_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void diskir_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xsdskb_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xsdili_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void diskm_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void disko_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void diskpbb_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xsdiskpn_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
// was  void xeq_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void C_equil(const double* energy, int nFlux, const double* params, int spectrumNumber, double* flux, double* fluxError, const char* initStr);
void xsxpdec_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void ezdiskbb_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xsgaul_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
// void xnneq_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void C_gnei(const double* energy, int nFlux, const double* params, int spectrumNumber, double* flux, double* fluxError, const char* initStr);
void grad_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xsgrbm_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void C_kerrbb(const double* energy, int nFlux, const double* params, int spectrumNumber, double* flux, double* fluxError, const char* initStr);
void C_kerrdisk(const double* energy, int nFlux, const double* params, int spectrumNumber, double* flux, double* fluxError, const char* initStr);
void spin_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void C_xslaor(const double* energy, int nFlux, const double* params, int spectrumNumber, double* flux, double* fluxError, const char* initStr);
// void laor2_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void C_laor2(const double* energy, int nFlux, const double* params, int spectrumNumber, double* flux, double* fluxError, const char* initStr);
void xslorz_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xsmeka_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xsmekl_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
// void xsmkcf_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void C_xsmkcf(const double* energy, int nFlux, const double* params, int spectrumNumber, double* flux, double* fluxError, const char* initStr);
// void C_xneq(const double* energy, int nFlux, const double* params, int spectrumNumber, double* flux, double* fluxError, const char* initStr);
void C_nei(const double* energy, int nFlux, const double* params, int spectrumNumber, double* flux, double* fluxError, const char* initStr);
// void xshock_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void C_npshock(const double* energy, int nFlux, const double* params, int spectrumNumber, double* flux, double* fluxError, const char* initStr);
void nsa_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void nsagrav_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void nsatmos_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void nsmax_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
// void xsnteea_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void C_xsnteea(const double* energy, int nFlux, const double* params, int spectrumNumber, double* flux, double* fluxError, const char* initStr);
// void nthcomp_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void C_nthcomp(const double* energy, int nFlux, const double* params, int spectrumNumber, double* flux, double* fluxError, const char* initStr);
void xspegp_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
// void xspexrav_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void C_xspexrav(const double* energy, int nFlux, const double* params, int spectrumNumber, double* flux, double* fluxError, const char* initStr);
void C_xspexriv(const double* energy, int nFlux, const double* params, int spectrumNumber, double* flux, double* fluxError, const char* initStr);
void xsp1tr_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void C_powerLaw(const double* energy, int nFlux, const double* params, int spectrumNumber, double* flux, double* fluxError, const char* initStr);
void xsposm_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
// void xneqs_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void C_pshock(const double* energy, int nFlux, const double* params, int spectrumNumber, double* flux, double* fluxError, const char* initStr);
void xsrays_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xredge_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xsrefsch_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
// void xsedov_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void C_sedov(const double* energy, int nFlux, const double* params, int spectrumNumber, double* flux, double* fluxError, const char* initStr);
void srcut_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void sresc_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xsstep_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xsvape_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xsbrmv_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
// void xseq_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void C_vequil(const double* energy, int nFlux, const double* params, int spectrumNumber, double* flux, double* fluxError, const char* initStr);
// void xsnneq_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void C_vgnei(const double* energy, int nFlux, const double* params, int spectrumNumber, double* flux, double* fluxError, const char* initStr);
void xsvmek_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xsvmkl_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
// void xsvmcf_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void C_xsvmcf(const double* energy, int nFlux, const double* params, int spectrumNumber, double* flux, double* fluxError, const char* initStr);
// void C_xsneq(const double* energy, int nFlux, const double* params, int spectrumNumber, double* flux, double* fluxError, const char* initStr);
void C_vnei(const double* energy, int nFlux, const double* params, int spectrumNumber, double* flux, double* fluxError, const char* initStr);
// void xsshock_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void C_vnpshock(const double* energy, int nFlux, const double* params, int spectrumNumber, double* flux, double* fluxError, const char* initStr);
// void xsneqs_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void C_vpshock(const double* energy, int nFlux, const double* params, int spectrumNumber, double* flux, double* fluxError, const char* initStr);
void xsvrys_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
// void xssedov_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void C_vsedov(const double* energy, int nFlux, const double* params, int spectrumNumber, double* flux, double* fluxError, const char* initStr);
void xszbod_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xszbrm_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
// void xszgau_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void C_xszgau(const double* energy, int nFlux, const double* params, int spectrumNumber, double* flux, double* fluxError, const char* initStr);
void C_zpowerLaw(const double* energy, int nFlux, const double* params, int spectrumNumber, double* flux, double* fluxError, const char* initStr);
void C_xsabsori(const double* energy, int nFlux, const double* params, int spectrumNumber, double* flux, double* fluxError, const char* initStr);
void acisabs_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xscnst_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xscabs_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xscycl_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xsdust_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xsedge_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xsabsc_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xsexp_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
// void xsgabs_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void C_gaussianAbsorptionLine(const double* energy, int nFlux, const double* params, int spectrumNumber, double* flux, double* fluxError, const char* initStr);
void xshecu_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xshrfl_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xsntch_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xsabsp_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xsphab_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xsplab_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
// void xspwab_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void C_xspwab(const double* energy, int nFlux, const double* params, int spectrumNumber, double* flux, double* fluxError, const char* initStr);
void xscred_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xssmdg_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void C_superExpCutoff(const double* energy, int nFlux, const double* params, int spectrumNumber, double* flux, double* fluxError, const char* initStr);
void xsspln_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xssssi_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void swind1_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
// void tbabs_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void C_tbabs(const double* energy, int nFlux, const double* params, int spectrumNumber, double* flux, double* fluxError, const char* initStr);
// void tbgrain_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void C_tbgrain(const double* energy, int nFlux, const double* params, int spectrumNumber, double* flux, double* fluxError, const char* initStr);
// void tbvabs_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void C_tbvabs(const double* energy, int nFlux, const double* params, int spectrumNumber, double* flux, double* fluxError, const char* initStr);
void xsred_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xsabsv_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xsvphb_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xsabsw_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xswnab_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xsxirf_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void mszdst_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xszedg_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xszhcu_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xszabp_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xszphb_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void zxipcf_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xszcrd_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void msldst_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
// void ztbabs_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void C_ztbabs(const double* energy, int nFlux, const double* params, int spectrumNumber, double* flux, double* fluxError, const char* initStr);
void xszvab_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xszvfe_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xszvph_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xszabs_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xszwnb_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);

// New XSPEC 12.7 models

void C_cplinear(const double* energy, int nFlux, const double* params, int spectrumNumber, double* flux, double* fluxError, const char* initStr);
// void xseqpair_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void C_xseqpair(const double* energy, int nFlux, const double* params, int spectrumNumber, double* flux, double* fluxError, const char* initStr);
// void xseqth_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void C_xseqth(const double* energy, int nFlux, const double* params, int spectrumNumber, double* flux, double* fluxError, const char* initStr);
// void xscompth_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void C_xscompth(const double* energy, int nFlux, const double* params, int spectrumNumber, double* flux, double* fluxError, const char* initStr);
void xsbvvp_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xsvvap_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void zigm_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);

// New XSPEC 12.7.1 models

void C_gaussDem(const double* energy, int nFlux, const double* params, int spectrumNumber, double* flux, double* fluxError, const char* initStr);
void C_vgaussDem(const double* energy, int nFlux, const double* params, int spectrumNumber, double* flux, double* fluxError, const char* initStr);
void logpar_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void eplogpar_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void optxagn_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void optxagnf_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void pexmon_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);

// Models from 12.8.0, 12.8.1, and 12.8.2

// additive
void C_agauss(const double* energy, int nFlux, const double* params, int spectrumNumber, double* flux, double* fluxError, const char* initStr);
void C_zagauss(const double* energy, int nFlux, const double* params, int spectrumNumber, double* flux, double* fluxError, const char* initStr);
void xscompmag(const double* energy, int nFlux, const double* params, int spectrumNumber, double* flux, double* fluxError, const char* initStr);
void xscomptb(const double* energy, int nFlux, const double* params, int spectrumNumber, double* flux, double* fluxError, const char* initStr);
void nsmaxg_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void nsx_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void C_rnei(const double* energy, int nFlux, const double* params, int spectrumNumber, double* flux, double* fluxError, const char* initStr);
void C_vrnei(const double* energy, int nFlux, const double* params, int spectrumNumber, double* flux, double* fluxError, const char* initStr);
void C_vvrnei(const double* energy, int nFlux, const double* params, int spectrumNumber, double* flux, double* fluxError, const char* initStr);
void C_vvgnei(const double* energy, int nFlux, const double* params, int spectrumNumber, double* flux, double* fluxError, const char* initStr);
void C_vvnei(const double* energy, int nFlux, const double* params, int spectrumNumber, double* flux, double* fluxError, const char* initStr);
void C_vvnpshock(const double* energy, int nFlux, const double* params, int spectrumNumber, double* flux, double* fluxError, const char* initStr);
void C_vvpshock(const double* energy, int nFlux, const double* params, int spectrumNumber, double* flux, double* fluxError, const char* initStr);
void C_vvsedov(const double* energy, int nFlux, const double* params, int spectrumNumber, double* flux, double* fluxError, const char* initStr);

//multiplicative
void xsphei_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xslyman_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xszbabs(const double* energy, int nFlux, const double* params, int spectrumNumber, double* flux, double* fluxError, const char* initStr);

// XSPEC table models
void xsatbl(float* ear, int ne, float* param, const char* filenm, int ifl, 
	    float* photar, float* photer);
void xsmtbl(float* ear, int ne, float* param, const char* filenm, int ifl, 
	    float* photar, float* photer);

// XSPEC convolution models
void C_cflux(const double* energy, int nFlux, const double* params, int spectrumNumber, double* flux, double* fluxError, const char* initStr);
void C_xsgsmt(const double* energy, int nFlux, const double* params, int spectrumNumber, double* flux, double* fluxError, const char* initStr);
void C_ireflct(const double* energy, int nFlux, const double* params, int spectrumNumber, double* flux, double* fluxError, const char* initStr);
void C_kdblur(const double* energy, int nFlux, const double* params, int spectrumNumber, double* flux, double* fluxError, const char* initStr);
void C_kdblur2(const double* energy, int nFlux, const double* params, int spectrumNumber, double* flux, double* fluxError, const char* initStr);
void C_spinconv(const double* energy, int nFlux, const double* params, int spectrumNumber, double* flux, double* fluxError, const char* initStr);
void C_xslsmt(const double* energy, int nFlux, const double* params, int spectrumNumber, double* flux, double* fluxError, const char* initStr);
void C_PartialCovering(const double* energy, int nFlux, const double* params, int spectrumNumber, double* flux, double* fluxError, const char* initStr);
void C_rdblur(const double* energy, int nFlux, const double* params, int spectrumNumber, double* flux, double* fluxError, const char* initStr);
void C_reflct(const double* energy, int nFlux, const double* params, int spectrumNumber, double* flux, double* fluxError, const char* initStr);
void C_simpl(const double* energy, int nFlux, const double* params, int spectrumNumber, double* flux, double* fluxError, const char* initStr);
void C_zashift(const double* energy, int nFlux, const double* params, int spectrumNumber, double* flux, double* fluxError, const char* initStr);
void C_zmshift(const double* energy, int nFlux, const double* params, int spectrumNumber, double* flux, double* fluxError, const char* initStr);
}

// Sun's C++ compiler complains if this is declared static
int _sherpa_init_xspec_library()
{

  static bool init = false;

  if ( init )
    return EXIT_SUCCESS;


  if ( !getenv("HEADAS") ) {
    // Raise appropriate error message that XSPEC initialization failed.
    PyErr_SetString( PyExc_ImportError,
		     (char*)"XSPEC initialization failed; "
		     "check HEADAS environment variable" );
    return EXIT_FAILURE;
  }
  
  // Stream buffer for redirected stdout
  std::streambuf* cout_sbuf = NULL;
  std::ofstream fout;

  try {

    // Redirect std::cout during call to FNINIT; this causes
    // us to swallow annoying messages about having set the
    // XSPEC abundance and cross section by default.

    // We have to do it this way because XSPEC code sends message
    // to XSPEC stream buffers, that in turn talk to std::cin,
    // std::cout, and std::cerr.  When we are not given a toggle
    // to turn off output to output streams, we can instead redirect
    // the stream itself.

    // In this case, the only stream we want to redirect is stdout.

    // If creating any stream buffer fails for any reason,
    // make sure we can still use original std::cout stream buffer.

    // Note:  since we are *not* redirecting std::cerr, any error
    // messages are still seen; if, for example, HEADAS is not set,
    // user will still see error message printed to stderr (i.e., to 
    // screen by default).  We still want such error messages to be
    // seen!  So do *not* redirect std::cerr.

    
    cout_sbuf = std::cout.rdbuf();
    fout.open("/dev/null");
    if (cout_sbuf != NULL && fout.is_open())
      std::cout.rdbuf(fout.rdbuf()); // temporary redirect stdout to /dev/null
    
    // Initialize XSPEC model library
    FNINIT();
    
    // Get back original std::cout
    if (cout_sbuf != NULL) {
      std::cout.clear();
      std::cout.rdbuf(cout_sbuf); 
    }
    fout.clear();
    fout.close();

    // Try to minimize model chatter for normal operation.
    FPCHAT( 0 );

    // Set cosmology initial values to XSPEC initial values
    csmph0( 70.0 );
    csmpq0( 0.0 );
    csmpl0( 0.73 );

  } catch(...) {

    
    // Get back original std::cout
    if (cout_sbuf != NULL) {
      std::cout.clear();
      std::cout.rdbuf(cout_sbuf);
    }
    fout.clear();
    fout.close();

    // Raise appropriate error message that XSPEC initialization failed.
    PyErr_SetString( PyExc_ImportError,
		     (char*)"XSPEC initialization failed; "
		     "check HEADAS environment variable" );
    return EXIT_FAILURE;

  }

  init = true;

  return EXIT_SUCCESS;

}

static PyObject* get_version( PyObject *self )
{ 

  if ( EXIT_SUCCESS != _sherpa_init_xspec_library() )
    return NULL;

  char version[256];
  int retval;

  try {
    
    retval = xs_getVersion(version, 256);

  } catch(...) {

    PyErr_SetString( PyExc_RuntimeError,
		     (char*)"could not get XSPEC version string" );
    return NULL;

  }
  
  if( retval < 0 ) {
    PyErr_SetString( PyExc_RuntimeError,
		     (char*)"XSPEC version string was truncated" );
    return NULL;
  }
  
  return Py_BuildValue( (char*)"s", version );

}

static PyObject* get_chatter( PyObject *self )
{ 

  if ( EXIT_SUCCESS != _sherpa_init_xspec_library() )
    return NULL;

  int chatter = 0;

  try {

    chatter = FGCHAT();

  } catch(...) {

    PyErr_SetString( PyExc_RuntimeError,
		     (char*)"could not get XSPEC chatter level" );
    return NULL;

  }

  return Py_BuildValue( (char*)"i", chatter );

}


static PyObject* get_abund( PyObject *self, PyObject *args )
{ 

  if ( EXIT_SUCCESS != _sherpa_init_xspec_library() )
    return NULL;

  char* abund = NULL;
  char* element = NULL;
  PyObject *retval = NULL;

  if ( !PyArg_ParseTuple( args, (char*)"|s", &element ) )
    return NULL;

  try {

    abund = FGSOLR();

  } catch(...) {

    PyErr_SetString( PyExc_RuntimeError,
		     (char*)"could not get XSPEC solar abundance" );
    return NULL;

  }

  if( !element ) {

    retval = (PyObject*) Py_BuildValue( (char*)"s", abund );
  
  } else {

    float abundVal = 0.0;
    std::streambuf *cerr_sbuf = NULL;
    std::ostringstream fout;
    
    try {
      
      cerr_sbuf = std::cerr.rdbuf();
      
      if (cerr_sbuf != NULL)
	std::cerr.rdbuf(fout.rdbuf());
      
      abundVal = FGABND(element);
      
      // Get back original std::cerr
      if (cerr_sbuf != NULL)
	std::cerr.rdbuf(cerr_sbuf); 
      
      
    } catch(...) {
      
      // Get back original std::cerr
      if (cerr_sbuf != NULL)
	std::cerr.rdbuf(cerr_sbuf);
      
      PyErr_Format( PyExc_RuntimeError,
		    (char*)"could not get XSPEC abundance for '%s'",
		  element);
      return NULL;
      
    }
    
    if( fout.str().size() > 0 ) {
      PyErr_Format( PyExc_TypeError,
		    (char*)"could not find element '%s'", element);
      return NULL;
    }
    
    retval = (PyObject*) Py_BuildValue( (char*)"f", abundVal );
  }

  return retval;
}


static PyObject* get_cosmo( PyObject *self )
{ 

  if ( EXIT_SUCCESS != _sherpa_init_xspec_library() )
    return NULL;

  float h0;
  float l0;
  float q0;

  try {

    h0 = csmgh0();
    l0 = csmgl0();
    q0 = csmgq0();

  } catch(...) {

    PyErr_SetString( PyExc_RuntimeError,
		     (char*)"could not get XSPEC cosmology settings" );
    return NULL;

  }

  return Py_BuildValue( (char*)"fff", h0, q0, l0 );

}


static PyObject* get_cross( PyObject *self )
{ 

  if ( EXIT_SUCCESS != _sherpa_init_xspec_library() )
    return NULL;

  char* cross = NULL;

  try {

    cross = FGXSCT();

  } catch(...) {

    PyErr_SetString( PyExc_RuntimeError,
		     (char*)"could not get XSPEC cross-section" );
    return NULL;

  }

  return Py_BuildValue( (char*)"s", cross );

}


static PyObject* set_chatter( PyObject *self, PyObject *args )
{ 

  if ( EXIT_SUCCESS != _sherpa_init_xspec_library() )
    return NULL;

  int chatter = 0;

  if ( !PyArg_ParseTuple( args, (char*)"i", &chatter ) )
    return NULL;

  try {

    FPCHAT( chatter );

  } catch(...) {
    
    PyErr_SetString( PyExc_RuntimeError,
		     (char*)"could not set XSPEC chatter level" );
    return NULL;

  }

  Py_RETURN_NONE;

}


static PyObject* set_abund( PyObject *self, PyObject *args )
{ 

  if ( EXIT_SUCCESS != _sherpa_init_xspec_library() )
    return NULL;

  char* table = NULL;
  int status = 0;

  if ( !PyArg_ParseTuple( args, (char*)"s", &table ) )
    return NULL;

  try {

    FPSOLR( table, &status );

  } catch(...) {
    
    status = 1;

  }

  // if abundance table name fails, try it as a filename
  if( status ) {

    std::ifstream fileStream(table);
    std::vector<float> vals(ABUND_SIZE, 0);
    size_t count(0);

    try {
      
      float element;
      fileStream.exceptions(std::ios_base::failbit);
      
      while (count < ABUND_SIZE && fileStream >> element) {
	vals[count] = element;
	++count;
      }
      
      status = 0;
    }
    catch ( std::exception& ) {
     
      if( !fileStream.eof() ) {
      	PyErr_Format( PyExc_ValueError, 
                      (char*)"Cannot read file '%s'.  It may not exist or contains invalid data",
                      table);
	return NULL;
      }
      
      status = 1;
    }

    try {
      
      FPSOLR((char*)"file", &status);
      FPSLFL( &vals[0], ABUND_SIZE, &status );
      
    } catch(...) {
      
      status = 1;
      
    }
  }

  if ( 0 != status ) {
    PyErr_SetString( PyExc_RuntimeError,
		     (char*)"could not set XSPEC abundance" );
    return NULL;
  }

  Py_RETURN_NONE;

}


static PyObject* set_cosmo( PyObject *self, PyObject *args )
{ 

  if ( EXIT_SUCCESS != _sherpa_init_xspec_library() )
    return NULL;

  float h0;
  float l0;
  float q0;

  if ( !PyArg_ParseTuple( args, (char*)"fff", &h0, &q0, &l0 ) )
    return NULL;

  try {

    csmph0( h0 );
    csmpl0( l0 );
    csmpq0( q0 );

  } catch(...) {

    PyErr_SetString( PyExc_RuntimeError,
		     (char*)"could not set XSPEC cosmology settings" );
    return NULL;

  }

  Py_RETURN_NONE;

}


static PyObject* set_cross( PyObject *self, PyObject *args )
{ 

  if ( EXIT_SUCCESS != _sherpa_init_xspec_library() )
    return NULL;

  char* csection = NULL;
  int status = 0;

  if ( !PyArg_ParseTuple( args, (char*)"s", &csection ) )
    return NULL;

  try {

    FPXSCT( csection, &status );

  } catch(...) {

    status = 1;

  }

  if ( 0 != status ) {
    PyErr_SetString( PyExc_RuntimeError,
		     (char*)"could not set XSPEC photoelectric "
		     "cross-section" );
    return NULL;
  }

  Py_RETURN_NONE;  

}


static PyObject* set_xset( PyObject *self, PyObject *args )
{ 

  if ( EXIT_SUCCESS != _sherpa_init_xspec_library() )
    return NULL;

  char* str_name = NULL;
  char* str_value = NULL;
  int status = 0;

  if ( !PyArg_ParseTuple( args, (char*)"ss", &str_name, &str_value ) )
    return NULL;

  try {

    FPMSTR( str_name, str_value );

  } catch(...) {

    status = 1;

  }

  if ( 0 != status ) {
    PyErr_Format( PyExc_RuntimeError,
		  (char*)"could not set XSPEC model strings '%s: %s'",
		  str_name, str_value);
    return NULL;
  }

  Py_RETURN_NONE;  

}

static PyObject* get_xset( PyObject *self, PyObject *args  )
{ 

  if ( EXIT_SUCCESS != _sherpa_init_xspec_library() )
    return NULL;

  char* str_name = NULL;
  char* str_value = NULL;

  if ( !PyArg_ParseTuple( args, (char*)"s", &str_name ) )
    return NULL;

  try {

    str_value = FGMSTR( str_name );

  } catch(...) {

    PyErr_Format( PyExc_RuntimeError,
		  (char*)"could not get XSPEC model string '%s'",
		  str_name);
    return NULL;

  }

  return Py_BuildValue( (char*)"s", str_value );

}

static PyMethodDef XSpecMethods[] = {
  { (char*)"get_xsversion", (PyCFunction)get_version, METH_NOARGS,
    (char*)"Get XSPEC version string." },

  { (char*)"get_xschatter", (PyCFunction)get_chatter, METH_NOARGS,
    (char*)"Get XSPEC chatter level." },
  { (char*)"set_xschatter", (PyCFunction)set_chatter, METH_VARARGS,
    (char*)"Set XSPEC chatter level." },

  { (char*)"get_xsabund", (PyCFunction)get_abund, METH_VARARGS,
    (char*)"Get XSPEC solar abundance." },
  { (char*)"set_xsabund", (PyCFunction)set_abund, METH_VARARGS,
    (char*)"Set XSPEC solar abundance." },

  { (char*)"set_xscosmo", (PyCFunction)set_cosmo, METH_VARARGS,
    (char*)"Set XSPEC cosmology settings (H0, q0, L0)." },
  { (char*)"get_xscosmo", (PyCFunction)get_cosmo, METH_NOARGS,
    (char*)"Get XSPEC cosmology settings (H0, q0, L0)." },

  { (char*)"get_xsxsect", (PyCFunction)get_cross, METH_NOARGS,
    (char*)"Get XSPEC photoelectric cross-section." },
  { (char*)"set_xsxsect", (PyCFunction)set_cross, METH_VARARGS,
    (char*)"Set XSPEC photoelectric cross-section." },

  { (char*)"set_xsxset", (PyCFunction)set_xset, METH_VARARGS,
    (char*)"Set XSPEC XSET <string name> <string value>" },
  { (char*)"get_xsxset", (PyCFunction)get_xset, METH_VARARGS,
    (char*)"Get XSPEC XSET <string name>" },

  XSPECMODELFCT_NORM( xsaped, 4 ),
  XSPECMODELFCT_NORM( xsbape, 5 ),
  XSPECMODELFCT_NORM( xsblbd, 2 ),
  XSPECMODELFCT_NORM( xsbbrd, 2 ),
  XSPECMODELFCT_C_NORM( C_xsbexrav, 10 ),
  XSPECMODELFCT_C_NORM( C_xsbexriv, 12 ),
  XSPECMODELFCT_C_NORM( C_brokenPowerLaw, 4 ),
  XSPECMODELFCT_C_NORM( C_broken2PowerLaw, 6 ),
  XSPECMODELFCT_C_NORM( C_sirf, 10 ),
  XSPECMODELFCT_NORM( xsbmc, 4 ),
  XSPECMODELFCT_NORM( xsbrms, 2 ),
  XSPECMODELFCT_NORM( xsbvpe, 17 ),
  XSPECMODELFCT_NORM( c6mekl, 11 ),
  XSPECMODELFCT_NORM( c6pmekl, 11 ),
  XSPECMODELFCT_NORM( c6pvmkl, 24 ),
  XSPECMODELFCT_NORM( c6vmekl, 24 ),
  XSPECMODELFCT_NORM( cemekl, 7 ),
  XSPECMODELFCT_C_NORM( C_cemVMekal, 20 ),
  XSPECMODELFCT_C_NORM( C_xscflw, 6 ),
  XSPECMODELFCT_NORM( compbb, 4 ),
  XSPECMODELFCT_NORM( compls, 3 ),
  XSPECMODELFCT_C_NORM( C_xscompps, 20 ),
  XSPECMODELFCT_NORM( compst, 3 ),
  XSPECMODELFCT_NORM( xstitg, 6 ),
  XSPECMODELFCT_C_NORM( C_cutoffPowerLaw, 3 ),
  XSPECMODELFCT_NORM( disk, 4 ),
  XSPECMODELFCT_NORM( diskir, 9 ),
  XSPECMODELFCT_NORM( xsdskb, 2 ),
  XSPECMODELFCT_NORM( xsdili, 6 ),
  XSPECMODELFCT_NORM( diskm, 5 ),
  XSPECMODELFCT_NORM( disko, 5 ),
  XSPECMODELFCT_NORM( diskpbb, 3 ),
  XSPECMODELFCT_NORM( xsdiskpn, 3 ),
  XSPECMODELFCT_C_NORM( C_equil, 4 ),
  XSPECMODELFCT_NORM( xsxpdec, 2 ),
  XSPECMODELFCT_NORM( ezdiskbb, 2 ),
  XSPECMODELFCT_NORM( xsgaul, 3 ),
  XSPECMODELFCT_C_NORM( C_gnei, 6 ),
  XSPECMODELFCT_NORM( grad, 7 ),
  XSPECMODELFCT_NORM( xsgrbm, 4 ),
  XSPECMODELFCT_C_NORM( C_kerrbb, 10 ),
  XSPECMODELFCT_C_NORM( C_kerrdisk, 8 ),
  XSPECMODELFCT_NORM( spin, 10 ),
  XSPECMODELFCT_C_NORM( C_xslaor, 6 ),
  XSPECMODELFCT_C_NORM( C_laor2, 8 ),
  XSPECMODELFCT_NORM( xslorz, 3 ),
  XSPECMODELFCT_NORM( xsmeka, 5 ),
  XSPECMODELFCT_NORM( xsmekl, 6 ),
  XSPECMODELFCT_C_NORM( C_xsmkcf, 6 ),
  XSPECMODELFCT_C_NORM( C_nei, 5 ),
  XSPECMODELFCT_C_NORM( C_npshock, 7 ),
  XSPECMODELFCT_NORM( nsa, 5 ),
  XSPECMODELFCT_NORM( nsagrav, 4 ),
  XSPECMODELFCT_NORM( nsatmos, 5 ),
  XSPECMODELFCT_NORM( nsmax, 4 ),
  XSPECMODELFCT_C_NORM( C_xsnteea, 16 ),
  XSPECMODELFCT_C_NORM( C_nthcomp, 6 ),
  XSPECMODELFCT_NORM( xspegp, 4 ),
  XSPECMODELFCT_C_NORM( C_xspexrav, 8 ),
  XSPECMODELFCT_C_NORM( C_xspexriv, 10 ),
  XSPECMODELFCT_NORM( xsp1tr, 11 ),
  XSPECMODELFCT_C_NORM( C_powerLaw, 2 ),
  XSPECMODELFCT_NORM( xsposm, 1 ),
  XSPECMODELFCT_C_NORM( C_pshock, 6 ),
  XSPECMODELFCT_NORM( xsrays, 4 ),
  XSPECMODELFCT_NORM( xredge, 3 ),
  XSPECMODELFCT_NORM( xsrefsch, 14 ),
  XSPECMODELFCT_C_NORM( C_sedov, 6 ),
  XSPECMODELFCT_NORM( srcut, 3 ),
  XSPECMODELFCT_NORM( sresc, 3 ),
  XSPECMODELFCT_NORM( xsstep, 3 ),
  XSPECMODELFCT_NORM( xsvape, 16 ),
  XSPECMODELFCT_NORM( xsbrmv, 3 ),
  XSPECMODELFCT_C_NORM( C_vequil, 15 ),
  XSPECMODELFCT_C_NORM( C_vgnei, 18 ),
  XSPECMODELFCT_NORM( xsvmek, 18 ),
  XSPECMODELFCT_NORM( xsvmkl, 19 ),
  XSPECMODELFCT_C_NORM( C_xsvmcf, 19 ),
  XSPECMODELFCT_C_NORM( C_vnei, 17 ),
  XSPECMODELFCT_C_NORM( C_vnpshock, 19 ),
  XSPECMODELFCT_C_NORM( C_vpshock, 18 ),
  XSPECMODELFCT_NORM( xsvrys, 15 ),
  XSPECMODELFCT_C_NORM( C_vsedov, 18 ),
  XSPECMODELFCT_NORM( xszbod, 3 ),
  XSPECMODELFCT_NORM( xszbrm, 3 ),
  XSPECMODELFCT_C_NORM( C_xszgau, 4 ),
  XSPECMODELFCT_C_NORM( C_zpowerLaw, 3 ),
  XSPECMODELFCT_C( C_xsabsori, 6 ),
  XSPECMODELFCT( acisabs, 8 ),
  XSPECMODELFCT( xscnst, 1 ),
  XSPECMODELFCT( xscabs, 1 ),
  XSPECMODELFCT( xscycl, 5 ),
  XSPECMODELFCT( xsdust, 2 ),
  XSPECMODELFCT( xsedge, 2 ),
  XSPECMODELFCT( xsabsc, 1 ),
  XSPECMODELFCT( xsexp, 3 ),
  XSPECMODELFCT_C( C_gaussianAbsorptionLine, 3 ),
  XSPECMODELFCT( xshecu, 2 ),
  XSPECMODELFCT( xshrfl, 8 ),
  XSPECMODELFCT( xsntch, 3 ),
  XSPECMODELFCT( xsabsp, 2 ),
  XSPECMODELFCT( xsphab, 1 ),
  XSPECMODELFCT( xsplab, 2 ),
  XSPECMODELFCT_C( C_xspwab, 3 ),
  XSPECMODELFCT( xscred, 1 ),
  XSPECMODELFCT( xssmdg, 4 ),
  XSPECMODELFCT_C( C_superExpCutoff, 2 ),
  XSPECMODELFCT( xsspln, 6 ),
  XSPECMODELFCT( xssssi, 1 ),
  XSPECMODELFCT( swind1, 4 ),
  XSPECMODELFCT_C( C_tbabs, 1 ),
  XSPECMODELFCT_C( C_tbgrain, 6 ),
  XSPECMODELFCT_C( C_tbvabs, 42 ),
  XSPECMODELFCT( xsred, 1 ),
  XSPECMODELFCT( xsabsv, 18 ),
  XSPECMODELFCT( xsvphb, 18 ),
  XSPECMODELFCT( xsabsw, 1 ),
  XSPECMODELFCT( xswnab, 2 ),
  XSPECMODELFCT( xsxirf, 13 ),
  XSPECMODELFCT( mszdst, 4 ),
  XSPECMODELFCT( xszedg, 3 ),
  XSPECMODELFCT( xszhcu, 3 ),
  XSPECMODELFCT( xszabp, 3 ),
  XSPECMODELFCT( xszphb, 2 ),
  XSPECMODELFCT( zxipcf, 4 ),
  XSPECMODELFCT( xszcrd, 2 ),
  XSPECMODELFCT( msldst, 4 ),
  XSPECMODELFCT_C( C_ztbabs, 2 ),
  XSPECMODELFCT( xszvab, 19 ),
  XSPECMODELFCT( xszvfe, 5 ),
  XSPECMODELFCT( xszvph, 19 ),
  XSPECMODELFCT( xszabs, 2 ),
  XSPECMODELFCT( xszwnb, 3 ),
  // New XSPEC 12.7 models
  XSPECMODELFCT_C_NORM( C_cplinear, 21 ),
  XSPECMODELFCT_C_NORM( C_xseqpair, 21 ),
  XSPECMODELFCT_C_NORM( C_xseqth, 21 ),
  XSPECMODELFCT_C_NORM( C_xscompth, 21 ),
  XSPECMODELFCT_NORM( xsbvvp, 34 ),
  XSPECMODELFCT_NORM( xsvvap, 33 ),
  XSPECMODELFCT( zigm, 3 ),
  // New XSPEC 12.7.1 models
  XSPECMODELFCT_C_NORM( C_gaussDem, 7 ),
  XSPECMODELFCT_C_NORM( C_vgaussDem, 20 ),
  XSPECMODELFCT_NORM( eplogpar, 3 ),
  XSPECMODELFCT_NORM( logpar, 4 ),
  XSPECMODELFCT_NORM( optxagn, 14 ),
  XSPECMODELFCT_NORM( optxagnf, 12 ),
  XSPECMODELFCT_NORM( pexmon, 8 ),

  // Models from 12.8.0, 12.8.1, and 12.8.2
  //
  // additive
  XSPECMODELFCT_C_NORM( C_agauss, 3 ),
  XSPECMODELFCT_C_NORM( C_zagauss, 4 ),
  XSPECMODELFCT_C_NORM( xscompmag, 9 ), // DJB thinks it's okay to use the C++ wrapper for C
  XSPECMODELFCT_C_NORM( xscomptb, 7 ), // DJB thinks it's okay to use the C++ wrapper for C
  XSPECMODELFCT_NORM( nsmaxg, 6 ),
  XSPECMODELFCT_NORM( nsx, 6 ),
  XSPECMODELFCT_C_NORM( C_rnei, 6 ),
  XSPECMODELFCT_C_NORM( C_vrnei, 18 ),
  XSPECMODELFCT_C_NORM( C_vvrnei, 35 ),
  XSPECMODELFCT_C_NORM( C_vvgnei, 35 ),
  XSPECMODELFCT_C_NORM( C_vvnei, 34 ),
  XSPECMODELFCT_C_NORM( C_vvnpshock, 36 ),
  XSPECMODELFCT_C_NORM( C_vvpshock, 35 ),
  XSPECMODELFCT_C_NORM( C_vvsedov, 35 ),

  //multiplicative
  XSPECMODELFCT( xsphei, 3 ),
  XSPECMODELFCT( xslyman, 4 ),
  XSPECMODELFCT_C( xszbabs, 4 ), // DJB thinks it's okay to use the C++ wrapper for C

  // XSPEC table models
  XSPECTABLEMODEL_NORM( xsatbl ),
  XSPECTABLEMODEL_NORM( xsmtbl ),
  // XSPEC convolution models
  XSPECMODELFCT_C(C_cflux, 3),
  XSPECMODELFCT_C(C_xsgsmt, 2),
  XSPECMODELFCT_C(C_ireflct, 7),
  XSPECMODELFCT_C(C_kdblur, 4),
  XSPECMODELFCT_C(C_kdblur2, 6),
  XSPECMODELFCT_C(C_spinconv, 7),
  XSPECMODELFCT_C(C_xslsmt, 2),
  XSPECMODELFCT_C(C_PartialCovering, 1),
  XSPECMODELFCT_C(C_rdblur, 4),
  XSPECMODELFCT_C(C_reflct, 5),
  XSPECMODELFCT_C(C_simpl, 3),
  XSPECMODELFCT_C(C_zashift, 1),
  XSPECMODELFCT_C(C_zmshift, 1),
  
  { NULL, NULL, 0, NULL }

};


PyMODINIT_FUNC
init_xspec(void)
{

  import_array();
  Py_InitModule( (char*)"_xspec", XSpecMethods );

}

