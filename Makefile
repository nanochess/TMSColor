# Ultra simple makefile for TMSColor
# by Oscar Toledo G.
# https://github.com/nanochess/tmscolor
#
build: tmscolor.o pletter.o
	@$(CC) tmscolor.o pletter.o -lm -o tmscolor

clean:
	@rm tmscolor.o pletter.o tmscolor

love:
	@echo "...not war"
