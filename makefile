## Main makefile.

include makefile_inc

# Rules:
# =====
default: mtklib tst exmples

all: mtklib tst exmples gendoc

exmples:
	@( cd $(EXAMPLES); $(MAKE) )

tst:
	@( cd $(TESTS); $(MAKE) )

checkheaders:
	@( cd $(INCLUDE); $(MAKE) )

mtklib:
	@( cd $(SRC); $(MAKE) )

gendoc:
	@( $(DOCGEN) $(DOCFILE) )
	@( cd $(DOC)/latex; $(MAKE) )

clean: cleanroot cleanlib cleanexamples cleanheaders

cleanroot:
	@(rm -f *~)

cleanlib:
	@( cd $(SRC); $(MAKE) clean )
	@( cd $(LIB); $(MAKE) clean )
	@( cd $(LIBDIR); rm -f libmtk.a )

cleantests:
	@( cd $(TESTS); $(MAKE) clean )

cleanexamples:
	@( cd $(EXAMPLES); $(MAKE) clean )

cleanheaders:
	@( cd $(INCLUDE); $(MAKE) clean)

release:
	@( $(MAKE) clean )
	@( cd ..; \
	   rm -f mtk*.tar.gz; \
	   tar czfv mtk_$(PLAT)-$(VERSION)-$(DATE).tar.gz mtk-v1.0/ )

bye: clean gendoc release

help:
	@echo '-----'
	@echo 'Makefile for MTK version 1.0.'
	@echo
	@echo 'Options are:'
	@echo '- make: builds only the library and the examples.'
	@echo '- all: builds the library, the examples and the documentation.'
	@echo '- mtklib: builds the library, i.e. generates the archive files.'
	@echo '- tests: generates the tests.'
	@echo '- examples: generates the examples.'
	@echo '- gendoc: generates the documentation for the library.'
	@echo '- checkheaders: checks syntax of the header files.'
	@echo
	@echo '- clean: cleans ALL the generated files.'
	@echo '- cleanlib: cleans the generated archive and object files.'
	@echo '- cleantests: cleans the generated tests executables.'
	@echo '- cleanexamples: cleans the generated examples executables.'
	@echo '-----'
