	' TMSColor 2.3 Jul/01/2024
	' Command: ./tmscolor -b -n -t1 statue3.bmp statue3.bas 
	' Created: Wed Sep 25 12:16:47 2024

	'
	' Recommended code:
	' MODE 0
	' DEFINE CHAR 1,13,image_char
	' DEFINE COLOR 1,13,image_color
	' SCREEN image_pattern,0,0,4,4,4
	'
	' Start tile = 1. Total_tiles = 13
image_char:
	DATA BYTE $ff,$ff,$ff,$ff,$ff,$ff,$ff,$ff
	DATA BYTE $ff,$ff,$ff,$ff,$f7,$f7,$ef,$ef
	DATA BYTE $ff,$ff,$ff,$df,$ef,$ef,$ef,$f7
	DATA BYTE $ff,$fe,$fc,$ff,$ff,$ff,$ff,$e0
	DATA BYTE $ef,$6f,$00,$48,$dd,$9f,$bf,$3f
	DATA BYTE $f7,$fb,$fb,$f9,$fc,$f8,$f8,$f8
	DATA BYTE $fd,$fa,$fe,$fe,$fe,$3e,$1e,$06
	DATA BYTE $00,$24,$7f,$7d,$18,$42,$e0,$f0
	DATA BYTE $ff,$ff,$ff,$ff,$ff,$ff,$7f,$1f
	DATA BYTE $f8,$fc,$fc,$fc,$fc,$fc,$fe,$ff
	DATA BYTE $00,$00,$00,$00,$00,$00,$00,$00
	DATA BYTE $f7,$e7,$0f,$00,$00,$01,$00,$00
	DATA BYTE $13,$27,$8f,$1f,$9f,$cf,$6f,$07

image_color:
	DATA BYTE $41,$41,$41,$41,$41,$41,$41,$41
	DATA BYTE $41,$41,$41,$41,$41,$41,$41,$41
	DATA BYTE $41,$41,$41,$41,$41,$41,$41,$41
	DATA BYTE $41,$41,$41,$41,$41,$41,$41,$41
	DATA BYTE $41,$41,$f1,$41,$41,$41,$41,$41
	DATA BYTE $41,$41,$41,$41,$41,$41,$41,$41
	DATA BYTE $41,$41,$41,$41,$41,$41,$41,$41
	DATA BYTE $f1,$41,$41,$41,$41,$41,$41,$41
	DATA BYTE $41,$41,$41,$41,$41,$41,$41,$41
	DATA BYTE $41,$41,$41,$41,$41,$41,$41,$41
	DATA BYTE $f1,$f1,$f1,$f1,$f1,$f1,$f1,$f1
	DATA BYTE $41,$41,$41,$f1,$f1,$41,$f1,$f1
	DATA BYTE $41,$41,$41,$41,$41,$41,$41,$41

	' Width = 4, height = 4
image_pattern:
	DATA BYTE $01,$01,$02,$01
	DATA BYTE $03,$04,$05,$01
	DATA BYTE $06,$07,$08,$09
	DATA BYTE $0a,$0b,$0c,$0d
