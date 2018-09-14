*****************************
The sherpa.astro.xspec module
*****************************

.. currentmodule:: sherpa.astro.xspec

.. automodule:: sherpa.astro.xspec

   .. rubric:: Classes
               
   .. autosummary::
      :toctree: api
   
      XSagauss
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
      XSbtapec
      XSbvapec
      XSbvtapec
      XSbvvapec
      XSbvvtapec
      XSc6mekl
      XSc6pmekl
      XSc6pvmkl
      XSc6vmekl
      XScarbatm
      XScemekl
      XScevmkl
      XScflow
      XScompLS
      XScompPS
      XScompST
      XScompTT
      XScompbb
      XScompmag
      XScomptb
      XScompth
      XScplinear
      XScutoffpl
      XSdisk
      XSdiskbb
      XSdiskir
      XSdiskline
      XSdiskm
      XSdisko
      XSdiskpbb
      XSdiskpn
      XSeplogpar
      XSeqpair
      XSeqtherm
      XSequil
      XSexpdec
      XSezdiskbb
      XSgadem
      XSgaussian
      XSgnei
      XSgrad
      XSgrbm
      XShatm
      XSkerrbb
      XSkerrd
      XSkerrdisk
      XSlaor
      XSlaor2
      XSlogpar
      XSlorentz
      XSmeka
      XSmekal
      XSmkcflow
      XSnei
      XSnlapec
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
      XSpegpwrlw
      XSpexmon
      XSpexrav
      XSpexriv
      XSplcabs
      XSposm
      XSpowerlaw
      XSpshock
      XSraymond
      XSredge
      XSrefsch
      XSrnei
      XSsedov
      XSsirf
      XSslimbh
      XSsnapec
      XSsrcut
      XSsresc
      XSstep
      XStapec
      XSvapec
      XSvbremss
      XSvequil
      XSvgadem
      XSvgnei
      XSvmcflow
      XSvmeka
      XSvmekal
      XSvnei
      XSvnpshock
      XSvoigt
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
      XSzagauss
      XSzbbody
      XSzbremss
      XSzgauss
      XSzpowerlw
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
      XScabs
      XSconstant
      XScyclabs
      XSdust
      XSedge
      XSexpabs
      XSexpfac
      XSgabs
      XSheilin
      XShighecut
      XShrefl
      XSismabs
      XSlyman
      XSnotch
      XSpcfabs
      XSphabs
      XSplabs
      XSpwab
      XSredden
      XSsmedge
      XSspexpcut
      XSspline
      XSswind1
      XSuvred
      XSvarabs
      XSvphabs
      XSwabs
      XSwndabs
      XSxion
      XSxscat
      XSzTBabs
      XSzbabs
      XSzdust
      XSzedge
      XSzhighect
      XSzigm
      XSzpcfabs
      XSzphabs
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

.. inheritance-diagram:: XSagauss XSapec XSbapec XSbbody XSbbodyrad XSbexrav XSbexriv XSbkn2pow XSbknpower XSbmc XSbremss XSbtapec XSbvapec XSbvtapec XSbvvapec XSbvvtapec XSc6mekl XSc6pmekl XSc6pvmkl XSc6vmekl XScarbatm XScemekl XScevmkl XScflow XScompLS XScompPS XScompST XScompTT XScompbb XScompmag XScomptb XScompth XScplinear XScutoffpl XSdisk XSdiskbb XSdiskir XSdiskline XSdiskm XSdisko XSdiskpbb XSdiskpn XSeplogpar XSeqpair XSeqtherm XSequil XSexpdec XSezdiskbb XSgadem XSgaussian XSgnei XSgrad XSgrbm XShatm XSkerrbb XSkerrd XSkerrdisk XSlaor XSlaor2 XSlogpar XSlorentz XSmeka XSmekal XSmkcflow XSnei XSnlapec XSnpshock XSnsa XSnsagrav XSnsatmos XSnsmax XSnsmaxg XSnsx XSnteea XSnthComp XSoptxagn XSoptxagnf XSpegpwrlw XSpexmon XSpexrav XSpexriv XSplcabs XSposm XSpowerlaw XSpshock XSraymond XSredge XSrefsch XSrnei XSsedov XSsirf XSslimbh XSsnapec XSsrcut XSsresc XSstep XStapec XSvapec XSvbremss XSvequil XSvgadem XSvgnei XSvmcflow XSvmeka XSvmekal XSvnei XSvnpshock XSvoigt XSvpshock XSvraymond XSvrnei XSvsedov XSvtapec XSvvapec XSvvgnei XSvvnei XSvvnpshock XSvvpshock XSvvrnei XSvvsedov XSvvtapec XSzagauss XSzbbody XSzbremss XSzgauss XSzpowerlw
   :parts: 1

.. inheritance-diagram:: XSSSS_ice XSTBabs XSTBfeo XSTBgas XSTBgrain XSTBpcf XSTBrel XSTBvarabs XSabsori XSacisabs XScabs XSconstant XScyclabs XSdust XSedge XSexpabs XSexpfac XSgabs XSheilin XShighecut XShrefl XSismabs XSlyman XSnotch XSpcfabs XSphabs XSplabs XSpwab XSredden XSsmedge XSspexpcut XSspline XSswind1 XSuvred XSvarabs XSvphabs XSwabs XSwndabs XSxion XSxscat XSzTBabs XSzbabs XSzdust XSzedge XSzhighect XSzigm XSzpcfabs XSzphabs XSzredden XSzsmdust XSzvarabs XSzvfeabs XSzvphabs XSzwabs XSzwndabs XSzxipcf
   :parts: 1


