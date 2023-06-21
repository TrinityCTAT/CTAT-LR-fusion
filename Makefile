all: ctat-minimap2_only plugins_only

ctat-minimap2_only:
	cd ctat-minimap2 && $(MAKE)

plugins_only:
	cd plugins && make

clean:
	cd ctat-minimap2 && $(MAKE) clean
	cd plugins && $(MAKE) clean

