SDDS_REPO = $(wildcard ../SDDS)
ifeq ($(SDDS_REPO),)
  $(error SDDS source code not found. Run 'git clone https://github.com/rtsoliday/SDDS.git' next to the SDDS-Python-module repository)
endif

include Makefile.rules

DIRS = $(SDDS_REPO)/zlib
DIRS += $(SDDS_REPO)/lzma
DIRS += $(SDDS_REPO)/mdblib
DIRS += $(SDDS_REPO)/mdbmth
DIRS += $(SDDS_REPO)/rpns/code
DIRS += $(SDDS_REPO)/namelist
DIRS += $(SDDS_REPO)/SDDSlib
DIRS += $(SDDS_REPO)/fftpack
DIRS += $(SDDS_REPO)/matlib
DIRS += $(SDDS_REPO)/mdbcommon
DIRS += src

.PHONY: all $(DIRS) clean distclean

all: $(DIRS)

$(SDDS_REPO)/zlib:
	$(MAKE) -C $@
$(SDDS_REPO)/lzma: $(SDDS_REPO)/zlib
	$(MAKE) -C $@
$(SDDS_REPO)/mdblib: $(SDDS_REPO)/lzma
	$(MAKE) -C $@
$(SDDS_REPO)/mdbmth: $(SDDS_REPO)/mdblib
	$(MAKE) -C $@
$(SDDS_REPO)/rpns/code: $(SDDS_REPO)/mdbmth
	$(MAKE) -C $@
$(SDDS_REPO)/namelist: $(SDDS_REPO)/rpns/code
	$(MAKE) -C $@
$(SDDS_REPO)/SDDSlib: $(SDDS_REPO)/namelist
	$(MAKE) -C $@
$(SDDS_REPO)/fftpack: $(SDDS_REPO)/SDDSlib
	$(MAKE) -C $@
$(SDDS_REPO)/matlib: $(SDDS_REPO)/fftpack
	$(MAKE) -C $@
$(SDDS_REPO)/mdbcommon: $(SDDS_REPO)/matlib
	$(MAKE) -C $@
$(SDDS_REPO)/utils: $(SDDS_REPO)/mdbcommon
	$(MAKE) -C $@
src: $(SDDS_REPO)/utils
	$(MAKE) -C $@

clean:
	$(MAKE) -C src clean

distclean: clean
	rm -rf bin/$(OS)-$(ARCH)
	rm -rf lib/$(OS)-$(ARCH)
