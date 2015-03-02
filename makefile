############################################################
FC = gfortran
SRCDIR = src
OBJDIR = obj
BINDIR = bin
INCDIR = inc
# vpath %.o $(OBJDIR)
# vpath %.f $(SRCDIR)
VPATH=$(SRCDIR)/main:$(OBJDIR):$(SRCDIR)/HelAmps:$(SRCDIR)/MatrixElements:\
      $(SRCDIR)/MultiCouplings:$(SRCDIR)/OffDiagonal:$(SRCDIR)/Required:\
      $(SRCDIR)/SMCorr:$(SRCDIR):$(SRCDIR)/ZZwidth
OUTPUT_OPTION= -o $(OBJDIR)/$@
FFLAGS = -g -Wall -ffixed-form -I$(INCDIR) -I$(INCDIR)/dists

############################################################
# main execs
MAIN = ffbarMultiZp-3G.o

# global dependencies

MADGRAPH = initialize.o alpha_EWNG.o 

PDFS = Cteq61Pdf.o setpdf.o

VEGAS = ve_dist.o rangen.o

HELAS = dhelas_all.o

# SMPROCESS = gg_ttb.o qqb_ttb.o bbb_ttb.o qqb_ttb_EW.o

SMPROCESS = gg_ttb.o qqb_ttb.o bbb_ttb.o 

FXN = fxn.o

FXNSYM = fxn_sym.o

# Z' widths: routine ZZWIDTHMULTI()

ZZWIDTH = ZZwidthMulti.o  

ZZWIDTHOD = ZZwidthMulti.o ZZwidthOD.o  

ZZWIDTHREAD = ZZwidth-Read.o  

ZZWIDTHODREAD = ZZwidth-Read.o ZZwidthOD-Read.o

ZZWIDTH4DCHM= ZZwidth-4DCHM.o

# Z' couplings: routine ZPCOUP()

ZPCOUP = MultiCouplings.o

ZPCOUP3G = MultiCouplings-3G.o

ZPCOUPAADDZA = MultiCouplings-AADD-ZA.o

ZPCOUPAADDA = MultiCouplings-AADD-A.o

ZPCOUPAADDZ = MultiCouplings-AADD-Z.o

ZPCOUPAADDWB = MultiCouplings-AADD-WB.o

ZPCOUPAADDB = MultiCouplings-AADD-B.o

ZPCOUPAADDW = MultiCouplings-AADD-W.o

ZPCOUPAADDOD = MultiCouplings-AADD-OD.o

ZPCOUPAADDSPLIT = MultiCouplings-AADD-Split.o

ZPCOUPAADDZAL1 = MultiCouplings-AADD-ZA-L1.o

ZPCOUPAADDWBL1 = MultiCouplings-AADD-WB-L1.o

ZPCOUPAADDFAKE = MultiCouplings-AADD-Fake.o 

ZPCOUPAADDFAKEWB = MultiCouplings-AADD-WB-Fake.o 

ZPCOUPZPNU = MultiCouplings-ZpNu.o ZZwidthOther.o

ZPCOUP4DCHM = MultiCouplings-4DCHM.o

# Corrections to Z boson fermionic couplings (e.g. due to Z-Z' mixing): routine SMCORR()

SMCORRNONE = SMcorr-none.o

SMCORR4DCHM = SMcorr-4DCHM.o

# f f~ > Z' > F F~ matrix element

BSMPROCESS =  qqb_ffb_multiZp_INT.o

BSMPROCESSOD =  QQ_2V_FF_OD.o diagonalise.o # off diagonal propagator version (2 Z' system only)

# determine q q~, g g > F F~ matrix elements

ZPME = HelAmps.o

ZPOD = HelAmps_OD.o

ZPONE = HelAmps_one.o # sets matrix element to 1 for checking phase space etc

############################################################
# Include prerequisites 
MAININC = vegas.inc ewparams.inc zpparams.inc runparams.inc names.inc VAcoups.inc \
LRcoups.inc dists_common.inc asy_common.inc plotparams.inc dists_bin.inc dists_plot.inc asy_plot.inc

FXNINC = vegas.inc ewparams.inc zpparams.inc runparams.inc dists_common.inc \
asy_common.inc asy_tot.inc dists_fill.inc asy_fill.inc

# Include files for plot management
DISTDIST = pT_dist.inc Et_dist.inc Pzff_dist.inc delY_dist.inc Yff_dist.inc beta_dist.inc \
Y3_dist.inc Y3col_dist.inc Y4_dist.inc Y4col_dist.inc \
cost3_dist.inc cost3col_dist.inc cost4_dist.inc cost4col_dist.inc \
eta3_dist.inc eta3col_dist.inc eta4_dist.inc eta4col_dist.inc

BINDIST = pT_bin.inc Et_bin.inc Pzff_bin.inc delY_bin.inc Yff_bin.inc beta_bin.inc \
Y3_bin.inc Y3col_bin.inc Y4_bin.inc Y4col_bin.inc \
cost3_bin.inc cost3col_bin.inc cost4_bin.inc cost4col_bin.inc \
eta3_bin.inc eta3col_bin.inc eta4_bin.inc eta4col_bin.inc

FILLDIST = pT_fill.inc Et_fill.inc Pzff_fill.inc delY_fill.inc Yff_fill.inc beta_fill.inc \
Y3_fill.inc Y3col_fill.inc Y4_fill.inc Y4col_fill.inc \
cost3_fill.inc cost3col_fill.inc cost4_fill.inc cost4col_fill.inc \
eta3_fill.inc eta3col_fill.inc eta4_fill.inc eta4col_fill.inc

PLOTDIST = pT_plot.inc Et_plot.inc Pzff_plot.inc delY_plot.inc Yff_plot.inc beta_plot.inc \
Y3_plot.inc Y3col_plot.inc Y4_plot.inc Y4col_plot.inc \
cost3_plot.inc cost3col_plot.inc cost4_plot.inc cost4col_plot.inc \
eta3_plot.inc eta3col_plot.inc eta4_plot.inc eta4col_plot.inc

FILLASY = ALL_fill.inc AL_fill.inc APV_fill.inc AFB_fill.inc AC_fill.inc \
AF_fill.inc AFB2_fill.inc AOFB_tot.inc ARFB_AOFB_fill.inc AFBSTAR_fill.inc

TOTASY = AFB_tot.inc AC_tot.inc AF_tot.inc AFB2_tot.inc ARFB_tot.inc AFBSTAR_tot.inc
############################################################
# Global dependencies
# GDEPS=$(MADGRAPH) $(VEGAS) $(HELAS) $(PDFS) $(SMPROCESS) $(patsubst %,inc/%, $(INCLUDE))
GDEPS=$(MADGRAPH) $(VEGAS) $(HELAS) $(PDFS) $(SMPROCESS) 
############################################################
EXEC = ffbarMultiZp ffbarRead ffbarODRead \
ffbarAADD_A ffbarAADD_Z ffbarAADD_B ffbarAADD_W ffbarAADD_WB ffbarAADD_ZA ffbarAADD_ZA_L1 \
ffbarAADD_WB_L1 ffbarAADD_Fake ffbarAADD_Split ffbarAADD_Split_OD ffbarAADD_OD ffbarAADD_Fake_WB \
ffbar4DCHM ffbarZpNu ffbarOne ffbarRead3G ffbarMultiZp3G
# HELOBJ = $(notdir, $($(wildcard src/HelAmps/*.f), .o=.f))

############################################################
.PHONY: clean
makefile: ;
############################################################
# manage include file dependencies
%.inc: 
	touch $@

$(INCDIR)/dists_common.inc: $(patsubst %,$(INCDIR)/dists/%, $(DISTDIST))

$(INCDIR)/dists_plot.inc: $(patsubst %,$(INCDIR)/dists/%, $(PLOTDIST))

$(INCDIR)/dists_fill.inc: $(patsubst %,$(INCDIR)/dists/%, $(FILLDIST))

$(INCDIR)/dists_bin.inc: $(patsubst %,$(INCDIR)/dists/%, $(BINDIST))

$(INCDIR)/asy_tot.inc: $(patsubst %,$(INCDIR)/dists/%, $(TOTASY))

$(INCDIR)/asy_fill.inc: $(patsubst %,$(INCDIR)/dists/%, $(FILLASY))

ffbarMultiZp-3G.o: $(patsubst %,$(INCDIR)/%, $(MAININC))

fxn.o fxn_sym.o: $(patsubst %,$(INCDIR)/%, $(FXNINC))

HELOBJ =$(notdir $(patsubst  %.f, %.o, $(wildcard $(SRCDIR)/HelAmps/*.f)))
$(HELOBJ): $(INCDIR)/runparams.inc $(INCDIR)/zpparams.inc $(INCDIR)/ewparams.inc

WIDOBJ =$(notdir $(patsubst  %.f, %.o, $(wildcard $(SRCDIR)/HZZwidth/*.f)))
$(WIDOBJ): $(INCDIR)/ewparams.inc $(INCDIR)/zpparams.inc $(INCDIR)/LRcoups.inc

ZZwidthRead.o ZZwidthODRead.o: names.inc

COUPOBJ =$(notdir $(patsubst  %.f, %.o, $(wildcard $(SRCDIR)/MultiCouplings/*.f)))
$(COUPOBJ): $(INCDIR)/VAcoups.inc

############################################################
all: $(patsubst %,$(BINDIR)/%, $(EXEC))

bin/%:
	$(eval DEPS=$(patsubst $(OBJDIR)/%,%,$^))
	$(F77) $(FFLAGS) -o $@ $(patsubst %,$(OBJDIR)/%, $(DEPS))
	
bin/ffbarMultiZp3G: $(MAIN) $(BSMPROCESS) $(ZZWIDTH) $(ZPCOUP3G) $(ZPME) $(SMCORRNONE) $(FXN) $(GDEPS)

bin/ffbarMultiZp3G_sym: $(MAIN) $(BSMPROCESS) $(ZZWIDTH) $(ZPCOUP3G) $(ZPME) $(SMCORRNONE) $(FXNSYM) $(GDEPS)

bin/ffbarMultiZpUniv: $(MAIN) $(BSMPROCESS) $(ZZWIDTH) $(ZPCOUP) $(ZPME) $(SMCORRNONE) $(FXN) $(GDEPS)

bin/ffbarOne: $(MAIN) $(BSMPROCESS) $(ZZWIDTH) $(ZPCOUP3G) $(ZPONE) $(SMCORRNONE) $(FXN) $(GDEPS)

bin/ffbarRead: $(MAIN) $(BSMPROCESS) $(ZZWIDTHREAD) $(ZPCOUP3G) $(ZPME) $(SMCORRNONE) $(FXN) $(GDEPS)

bin/ffbarODRead: $(MAIN) $(BSMPROCESSOD) $(ZZWIDTHODREAD) $(ZPCOUP3G) $(ZPOD) $(SMCORRNONE) $(FXN) $(GDEPS)

bin/ffbar4DCHM: $(MAIN) $(BSMPROCESS) $(ZPCOUP4DCHM) $(ZZWIDTH4DCHM) $(ZPME) $(SMCORR4DCHM) $(FXN) $(GDEPS)

bin/ffbarZpNu: $(MAIN) $(BSMPROCESS) $(ZZWIDTH) $(ZPCOUPZPNU) $(ZPME) $(SMCORRNONE) $(FXN) $(GDEPS)

bin/ffbarAADD_ZA: $(MAIN) $(BSMPROCESS) $(ZZWIDTH) $(ZPCOUPAADDZA) $(ZPME) $(SMCORRNONE) $(FXN) $(GDEPS)

bin/ffbarAADD_WB: $(MAIN) $(BSMPROCESS) $(ZZWIDTH) $(ZPCOUPAADDWB) $(ZPME) $(SMCORRNONE) $(FXN) $(GDEPS)

bin/ffbarAADD_B: $(MAIN) $(BSMPROCESS) $(ZZWIDTH) $(ZPCOUPAADDB) $(ZPME) $(SMCORRNONE) $(FXN) $(GDEPS)

bin/ffbarAADD_W: $(MAIN) $(BSMPROCESS) $(ZZWIDTH) $(ZPCOUPAADDW) $(ZPME) $(SMCORRNONE) $(FXN) $(GDEPS)

bin/ffbarAADD_A: $(MAIN) $(BSMPROCESS) $(ZZWIDTH) $(ZPCOUPAADDA) $(ZPME) $(SMCORRNONE) $(FXN) $(GDEPS)

bin/ffbarAADD_Z: $(MAIN) $(BSMPROCESS) $(ZZWIDTH) $(ZPCOUPAADDZ) $(ZPME) $(SMCORRNONE) $(FXN) $(GDEPS)

bin/ffbarAADD_ZA_L1: $(MAIN) $(BSMPROCESS) $(ZZWIDTH) $(ZPCOUPAADDZAL1) $(ZPME) $(SMCORRNONE) $(FXN) $(GDEPS)

bin/ffbarAADD_WB_L1: $(MAIN) $(BSMPROCESS) $(ZZWIDTH) $(ZPCOUPAADDWBL1) $(ZPME) $(SMCORRNONE) $(FXN) $(GDEPS)

bin/ffbarAADD_Fake: $(MAIN) $(BSMPROCESS) $(ZPCOUPAADDFAKE) $(ZPME) $(SMCORRNONE) $(FXN) $(GDEPS)

bin/ffbarAADD_Fake_WB: $(MAIN) $(BSMPROCESS) $(ZPCOUPAADDFAKEWB) $(ZPME) $(SMCORRNONE) $(FXN) $(GDEPS)

bin/ffbarAADD_OD: $(MAIN) $(BSMPROCESSOD) $(ZZWIDTHOD) $(ZPCOUPAADDZA) $(ZPOD) $(SMCORRNONE) $(FXN) $(GDEPS)

bin/ffbarAADD_Split: $(MAIN) $(BSMPROCESS) $(ZZWIDTH) $(ZPCOUPAADDSPLIT) $(ZPME) $(SMCORRNONE) $(FXN) $(GDEPS)

bin/ffbarAADD_Split_OD: $(MAIN) $(BSMPROCESSOD) $(ZZWIDTHOD) $(ZPCOUPAADDSPLIT) $(ZPOD) $(SMCORRNONE) $(FXN) $(GDEPS)

clean:
	rm -f $(OBJDIR)/*.o 
	rm -f $(BINDIR)/* 

