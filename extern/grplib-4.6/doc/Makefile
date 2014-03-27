#-----------------------------------------------------------------------
# Makefile for the pipe/tools documentation.
#
# To make a postscript version of the documentation type:
#       % make
#
# To install the documentation type:
#       % make install
#
#-----------------------------------------------------------------------
MK_TOP = ../../..
include $(MK_TOP)/Makefile.master
include $(MK_TOP)/include/Makefile.scidev

PORT_XML	  = \
	group.xml \
	grpadaptive.xml \
	grpadaptivesnr.xml \
	grpbin.xml \
	grpbinfile.xml \
	grpbinwidth.xml \
	grpgetchanspergroup.xml \
	grpgetgroupsum.xml \
	grpgetgrpnum.xml \
	grpmaxslope.xml \
	grpminslope.xml \
	grpnumbins.xml \
	grpnumcounts.xml \
	grpsnr.xml 

include $(MK_TOP)/Makefile.all

#-----------------------------------------------------------------------
#                       MAKEFILE DEPENDENCIES
#-----------------------------------------------------------------------

announce1:
	@echo "   /----------------------------------------------------------\ "
	@echo "   |              Building pyglue documentation               | "
	@echo "   \----------------------------------------------------------/ "
