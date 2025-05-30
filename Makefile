# Ultra simple makefile for TMSColor
# by Oscar Toledo G.
# https://github.com/nanochess/tmscolor
#
tmscolor: tmscolor.o pletter.o lodepng.o
	@$(CC) tmscolor.o pletter.o lodepng.o -lm -o $@

clean:
	@rm tmscolor.o pletter.o lodepng.o tmscolor

love:
	@echo "...not war"
