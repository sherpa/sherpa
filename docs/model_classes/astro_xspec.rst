*****************************
The sherpa.astro.xspec module
*****************************

.. currentmodule:: sherpa.astro.xspec

.. automodule:: sherpa.astro.xspec

   .. rubric:: Classes

   .. autosummary::
      :toctree: api

      XSSSS_ice
      XSTBabs
      XSTBfeo
      XSTBgas
      XSTBgrain
      XSTBpcf
      XSTBrel
      XSTBvarabs
      XSabsori
      XSacisabs
      XSagauss
      XSagnsed
      XSapec
      XSbapec
      XSbbody
      XSbbodyrad
      XSbexrav
      XSbexriv
      XSbkn2pow
      XSbknpower
      XSbmc
      XSbremss
      XSbrnei
      XSbtapec
      XSbvapec
      XSbvrnei
      XSbvtapec
      XSbvvapec
      XSbvvrnei
      XSbvvtapec
      XSbwcycl
      XSc6mekl
      XSc6pmekl
      XSc6pvmkl
      XSc6vmekl
      XScabs
      XScarbatm
      XScemekl
      XScevmkl
      XScflow
      XScflux
      XScglumin
      XSclumin
      XScompLS
      XScompPS
      XScompST
      XScompTT
      XScompbb
      XScompmag
      XScomptb
      XScompth
      XSconstant
      XScpflux
      XScph
      XScplinear
      XScutoffpl
      XScyclabs
      XSdisk
      XSdiskbb
      XSdiskir
      XSdiskline
      XSdiskm
      XSdisko
      XSdiskpbb
      XSdiskpn
      XSdust
      XSedge
      XSeplogpar
      XSeqpair
      XSeqtherm
      XSequil
      XSexpabs
      XSexpdec
      XSexpfac
      XSezdiskbb
      XSgabs
      XSgadem
      XSgaussian
      XSgnei
      XSgrad
      XSgrbcomp
      XSgrbm
      XSgsmooth
      XShatm
      XSheilin
      XShighecut
      XShrefl
      XSireflect
      XSismabs
      XSjet
      XSkdblur
      XSkdblur2
      XSkerrbb
      XSkerrconv
      XSkerrd
      XSkerrdisk
      XSkyconv
      XSkyrline
      XSlaor
      XSlaor2
      XSlogpar
      XSlorentz
      XSlsmooth
      XSlyman
      XSmeka
      XSmekal
      XSmkcflow
      XSnei
      XSnlapec
      XSnotch
      XSnpshock
      XSnsa
      XSnsagrav
      XSnsatmos
      XSnsmax
      XSnsmaxg
      XSnsx
      XSnteea
      XSnthComp
      XSoptxagn
      XSoptxagnf
      XSpartcov
      XSpcfabs
      XSpegpwrlw
      XSpexmon
      XSpexrav
      XSpexriv
      XSphabs
      XSplabs
      XSplcabs
      XSposm
      XSpowerlaw
      XSpshock
      XSpwab
      XSqsosed
      XSraymond
      XSrdblur
      XSredden
      XSredge
      XSreflect
      XSrefsch
      XSrfxconv
      XSrgsxsrc
      XSrnei
      XSsedov
      XSsimpl
      XSsirf
      XSslimbh
      XSsmedge
      XSsnapec
      XSspexpcut
      XSspline
      XSsrcut
      XSsresc
      XSssa
      XSstep
      XSswind1
      XStapec
      XSthcomp
      XSuvred
      XSvapec
      XSvarabs
      XSvashift
      XSvbremss
      XSvcph
      XSvequil
      XSvgadem
      XSvgnei
      XSvmcflow
      XSvmeka
      XSvmekal
      XSvmshift
      XSvnei
      XSvnpshock
      XSvoigt
      XSvphabs
      XSvpshock
      XSvraymond
      XSvrnei
      XSvsedov
      XSvtapec
      XSvvapec
      XSvvgnei
      XSvvnei
      XSvvnpshock
      XSvvpshock
      XSvvrnei
      XSvvsedov
      XSvvtapec
      XSwabs
      XSwndabs
      XSxilconv
      XSxion
      XSxscat
      XSzTBabs
      XSzagauss
      XSzashift
      XSzbabs
      XSzbbody
      XSzbknpower
      XSzbremss
      XSzcutoffpl
      XSzdust
      XSzedge
      XSzgauss
      XSzhighect
      XSzigm
      XSzlogpar
      XSzmshift
      XSzpcfabs
      XSzphabs
      XSzpowerlw
      XSzredden
      XSzsmdust
      XSzvarabs
      XSzvfeabs
      XSzvphabs
      XSzwabs
      XSzwndabs
      XSzxipcf

Class Inheritance Diagram
=========================

.. inheritance-diagram:: XSagauss XSagnsed XSagnslim XSapec XSbapec XSbbody XSbbodyrad XSbexrav XSbexriv XSbkn2pow XSbknpower XSbmc XSbremss XSbrnei XSbtapec XSbvapec XSbvrnei XSbvtapec XSbvvapec XSbvvrnei XSbvvtapec XSbwcycl XSc6mekl XSc6pmekl XSc6pvmkl XSc6vmekl XScarbatm XScemekl XScevmkl XScflow XScompLS XScompPS XScompST XScompTT XScompbb XScompmag XScomptb XScompth XScph XScplinear XScutoffpl XSdisk XSdiskbb XSdiskir XSdiskline XSdiskm XSdisko XSdiskpbb XSdiskpn XSeplogpar XSeqpair XSeqtherm XSequil XSexpdec XSezdiskbb XSgadem XSgaussian XSgnei XSgrad XSgrbcomp XSgrbm XShatm XSjet XSkerrbb XSkerrd XSkerrdisk XSkyrline XSlaor XSlaor2 XSlogpar XSlorentz XSmeka XSmekal XSmkcflow XSnei XSnlapec XSnpshock XSnsa XSnsagrav XSnsatmos XSnsmax XSnsmaxg XSnsx XSnteea XSnthComp XSoptxagn XSoptxagnf XSpegpwrlw XSpexmon XSpexrav XSpexriv XSplcabs XSposm XSpowerlaw XSpshock XSqsosed XSraymond XSredge XSrefsch XSrnei XSsedov XSsirf XSslimbh XSsnapec XSsrcut XSsresc XSssa XSstep XStapec XSvapec XSvbremss XSvcph XSvequil XSvgadem XSvgnei XSvmcflow XSvmeka XSvmekal XSvnei XSvnpshock XSvoigt XSvpshock XSvraymond XSvrnei XSvsedov XSvtapec XSvvapec XSvvgnei XSvvnei XSvvnpshock XSvvpshock XSvvrnei XSvvsedov XSvvtapec XSzagauss XSzbbody XSzbknpower XSzbremss XSzcutoffpl XSzgauss XSzkerrbb XSzlogpar XSzpowerlw
   :parts: 1

.. inheritance-diagram:: XSSSS_ice XSTBabs XSTBfeo XSTBgas XSTBgrain XSTBpcf XSTBrel XSTBvarabs XSabsori XSacisabs XScabs XSconstant XScyclabs XSdust XSedge XSexpabs XSexpfac XSgabs XSheilin XShighecut XShrefl XSismabs XSismdust XSlog10con XSlogconst XSlyman XSnotch XSolivineabs XSpcfabs XSphabs XSplabs XSpwab XSredden XSsmedge XSspexpcut XSspline XSswind1 XSuvred XSvarabs XSvphabs XSwabs XSwndabs XSxion XSxscat XSzTBabs XSzbabs XSzdust XSzedge XSzhighect XSzigm XSzpcfabs XSzphabs XSzredden XSzsmdust XSzvarabs XSzvfeabs XSzvphabs XSzwabs XSzwndabs XSzxipcf
   :parts: 1

.. inheritance-diagram:: XScflux XScglumin XSclumin XScpflux XSgsmooth XSireflect XSkdblur XSkdblur2 XSkerrconv XSkyconv XSlsmooth XSpartcov XSrdblur XSreflect XSrfxconv XSrgsxsrc XSsimpl XSthcomp XSvashift XSvmshift XSxilconv XSzashift XSzmshift
   :parts: 1
