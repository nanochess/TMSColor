# Ultra simple makefile for TMSColor
# by Oscar Toledo G.
# https://github.com/nanochess/tmscolor
#
build:
	@cc tmscolor.c -lm -o tmscolor

clean:
	@rm tmscolor

love:
	@echo "...not war"

