/*
 ** TMSColor: Converts a BMP image to TMS9928 bitmap/color format
 **
 ** by Oscar Toledo Gutiérrez
 ** http://nanochess.org/
 **
 ** Copyright (C) 2009-2025 Oscar Toledo Gutiérrez
 **
 ** This program is free software; you can redistribute it and/or modify
 ** it under the terms of the GNU General Public License as published by
 ** the Free Software Foundation; either version 2 of the License, or
 ** (at your option) any later version.
 **
 ** This program is distributed in the hope that it will be useful,
 ** but WITHOUT ANY WARRANTY; without even the implied warranty of
 ** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 ** GNU General Public License for more details.
 **
 ** You should have received a copy of the GNU General Public License along
 ** with this program; if not, write to the Free Software Foundation, Inc.,
 ** 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 **
 **
 ** Creation date: Jun/06/2009.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include <limits.h>
#include <math.h>
#include "lodepng.h"

#define VERSION "3.2 Aug/26/2025"     /* Software version */

#define ROUND8(x)  ((x + 7) & ~7)

extern void pletter(unsigned char *, int, unsigned char **, int *);

unsigned char *image;           /* Image data in RGB format */
int size_x;                     /* Size X in pixels */
int size_y;                     /* Size Y in pixels */

unsigned char *bitmap;
unsigned char *color;
unsigned char *pattern;
unsigned char *source;

unsigned char *source2;
unsigned char *ignore;
static unsigned char usage[24][32][32];
unsigned char sprites[2048];
unsigned char attr[128];

int use_bitmap = 0;
int sms_mode = 0;
int nes_mode = 0;

/*
 ** Use this palette in your paint program
 */
unsigned char colors[16 * 3] = {
    0x40, 0x40, 0x40,
    0x00, 0x00, 0x00,
    0x45, 0xc5, 0x25,    /* 2 - Green */
    0x7a, 0xda, 0x62,    /* 3 - Light green */
    0xe8, 0x57, 0x55,    /* 4 - Blue */
    0xf8, 0x77, 0x7d,    /* 5 - Light blue */
    0x4e, 0x53, 0xd0,    /* 6 - Dark red */
    0xf2, 0xea, 0x47,    /* 7 - Cyan */
    0x56, 0x57, 0xf8,    /* 8 - Red */
    0x7a, 0x7b, 0xff,    /* 9 - Light red */
    0x58, 0xc0, 0xd3,    /* 10 - Yellow */
    0x83, 0xcd, 0xe5,    /* 11 - Light yellow */
    0x3d, 0xad, 0x25,    /* 12 - Dark green */
    0xb8, 0x5d, 0xc7,
    0xcc, 0xcc, 0xcc,    /* 14 - Gray */
    0xff, 0xff, 0xff,    /* 15 - White */
};

int flip_bits[256];

/*
 ** Converts from hexadecimal
 */
int from_hex(int letter)
{
    letter = toupper(letter);
    if (letter < '0')
        return 0;
    if (letter > 'F')
        return 15;
    letter -= '0';
    if (letter > 9)
        letter -= 7;
    return letter;
}

/*
 ** Prototypes
 */
int main(int, char *[]);

int check_flip_x(unsigned char *saved, unsigned char *new)
{
    int c;
    
    for (c = 0; c < 32; c++) {
        if (flip_bits[saved[c]] != new[c])
            return 0;
    }
    return 1;
}

int check_flip_y(unsigned char *saved, unsigned char *new)
{
    int c;
    
    for (c = 0; c < 32; c += 4) {
        if (saved[28 - c] != new[c])
            return 0;
        if (saved[29 - c] != new[c + 1])
            return 0;
        if (saved[30 - c] != new[c + 2])
            return 0;
        if (saved[31 - c] != new[c + 3])
            return 0;
    }
    return 1;
}

int check_flip_xy(unsigned char *saved, unsigned char *new)
{
    int c;
    
    for (c = 0; c < 32; c += 4) {
        if (flip_bits[saved[28 - c]] != new[c])
            return 0;
        if (flip_bits[saved[29 - c]] != new[c + 1])
            return 0;
        if (flip_bits[saved[30 - c]] != new[c + 2])
            return 0;
        if (flip_bits[saved[31 - c]] != new[c + 3])
            return 0;
    }
    return 1;
}

int check_triple_color(void)
{
    int c;
    int d;
    int y;
    int x;
    int g;
    int b;
    int color2;
    int color3;
    int triple_color;
    
    memset(&usage[0][0][0], 255, sizeof(usage));
    memset(ignore, 0, size_x * size_y);
    //
    // Find the color usage per 8x8 block
    //
    triple_color = 0;
    for (c = 0; c < 24; c++) {
        for (d = 0; d < 32; d++) {
            for (y = c * 8; y < (c + 1) * 8; y++) {
                color2 = 255;
                color3 = 255;
                g = 255;
                for (x = d * 8; x < (d + 1) * 8; x++) {
                    if (color2 == source[y * size_x + x])
                        ;
                    else if (color3 == source[y * size_x + x])
                        ;
                    else if (color2 == 255)
                        color2 = source[y * size_x + x];
                    else if (color3 == 255)
                        color3 = source[y * size_x + x];
                    else {
                        if (g == 255)
                            triple_color++;
                        g = source[y * size_x + x];
                        for (b = 0; b < 16; b++) {
                            if (usage[c][d][b] == color2)
                                break;
                            if (usage[c][d][b] == 255)
                                break;
                        }
                        if (usage[c][d][b] == 255) {
                            usage[c][d][b] = color2;
                            usage[c][d][color2 + 16] = 0;
                        }
                        for (b = 0; b < 16; b++) {
                            if (usage[c][d][b] == color3)
                                break;
                            if (usage[c][d][b] == 255)
                                break;
                        }
                        if (usage[c][d][b] == 255) {
                            usage[c][d][b] = color3;
                            usage[c][d][color3 + 16] = 0;
                        }
                        for (b = 0; b < 16; b++) {
                            if (usage[c][d][b] == source[y * size_x + x])
                                break;
                            if (usage[c][d][b] == 255)
                                break;
                        }
                        if (usage[c][d][b] == 255) {
                            usage[c][d][b] = source[y * size_x + x];
                            usage[c][d][source[y * size_x + x] + 16] = 0;
                        }
                    }
                }
                if (g == 255) {
                    for (x = d * 8; x < (d + 1) * 8; x++)
                        ignore[y * size_x + x] = 1;
                }
            }
        }
    }
    return triple_color;
}

#define sqr(c)  ((c) * (c))

typedef struct {
    double a;
    double b;
    double l;
} LAB;

/*
 ** Converted from Javascript code at https://github.com/antimatter15/rgb-lab
 ** MIT license Copyright (c) 2014 Kevin Kwok <antimatter15@gmail.com>
 */
void rgb2lab(unsigned char rgb[3], LAB *lab)
{
    double r = rgb[2] / 255.0;
    double g = rgb[1] / 255.0;
    double b = rgb[0] / 255.0;
    double x;
    double y;
    double z;
    
    r = (r > 0.04045) ? pow((r + 0.055) / 1.055, 2.4) : r / 12.92;
    g = (g > 0.04045) ? pow((g + 0.055) / 1.055, 2.4) : g / 12.92;
    b = (b > 0.04045) ? pow((b + 0.055) / 1.055, 2.4) : b / 12.92;
    
    x = (r * 0.4124 + g * 0.3576 + b * 0.1805) / 0.95047;
    y = (r * 0.2126 + g * 0.7152 + b * 0.0722) / 1.00000;
    z = (r * 0.0193 + g * 0.1192 + b * 0.9505) / 1.08883;
    
    x = (x > 0.008856) ? pow(x, 1.0/3) : (7.787 * x) + 16.0/116;
    y = (y > 0.008856) ? pow(y, 1.0/3) : (7.787 * y) + 16.0/116;
    z = (z > 0.008856) ? pow(z, 1.0/3) : (7.787 * z) + 16.0/116;
    
    lab->l = (116 * y) - 16;
    lab->a = 500 * (x - y);
    lab->b = 200 * (y - z);
}

/*
 ** Converted from C++ code at https://github.com/gfiumara/CIEDE2000
 ** MIT license Copyright (c) 2015 Greg Fiumara
 */
#ifndef M_PI
#define M_PI 3.14159265358979
#endif

#define deg2Rad(deg)    ((deg) * (M_PI / 180.0))

#define rad2Deg(rad)    ((180.0 / M_PI) * (rad))

#define k_L 1.0
#define k_C 1.0
#define k_H 1.0
#define deg360InRad deg2Rad(360.0)
#define deg180InRad deg2Rad(180.0)
#define pow25To7    6103515625.0    /* pow(25, 7) */

double CIEDE2000(LAB *lab1, LAB *lab2)
{
    
    /*
     * Step 1
     */
    /* Equation 2 */
    double C1 = sqrt((lab1->a * lab1->a) + (lab1->b * lab1->b));
    double C2 = sqrt((lab2->a * lab2->a) + (lab2->b * lab2->b));
    /* Equation 3 */
    double barC = (C1 + C2) / 2.0;
    /* Equation 4 */
    double G = 0.5 * (1 - sqrt(pow(barC, 7) / (pow(barC, 7) + pow25To7)));
    /* Equation 5 */
    double a1Prime = (1.0 + G) * lab1->a;
    double a2Prime = (1.0 + G) * lab2->a;
    /* Equation 6 */
    double CPrime1 = sqrt((a1Prime * a1Prime) + (lab1->b * lab1->b));
    double CPrime2 = sqrt((a2Prime * a2Prime) + (lab2->b * lab2->b));
    /* Equation 7 */
    double hPrime1;
    double hPrime2;
    double deltaLPrime;
    double deltaCPrime;
    double deltahPrime;
    double CPrimeProduct;
    double deltaHPrime;
    double barLPrime;
    double barCPrime;
    double barhPrime;
    double hPrimeSum;
    double T;
    double deltaTheta;
    double R_C;
    double S_L;
    double S_C;
    double S_H;
    double R_T;
    double deltaE;
    
    if (lab1->b == 0 && a1Prime == 0)
        hPrime1 = 0.0;
    else {
        hPrime1 = atan2(lab1->b, a1Prime);
        if (hPrime1 < 0)
            hPrime1 += deg360InRad;
    }
    if (lab2->b == 0 && a2Prime == 0)
        hPrime2 = 0.0;
    else {
        hPrime2 = atan2(lab2->b, a2Prime);
        if (hPrime2 < 0)
            hPrime2 += deg360InRad;
    }
    
    /*
     * Step 2
     */
    /* Equation 8 */
    deltaLPrime = lab2->l - lab1->l;
    /* Equation 9 */
    deltaCPrime = CPrime2 - CPrime1;
    /* Equation 10 */
    
    CPrimeProduct = CPrime1 * CPrime2;
    if (CPrimeProduct == 0)
        deltahPrime = 0;
    else {
        /* Avoid the fabs() call */
        deltahPrime = hPrime2 - hPrime1;
        if (deltahPrime < -deg180InRad)
            deltahPrime += deg360InRad;
        else if (deltahPrime > deg180InRad)
            deltahPrime -= deg360InRad;
    }
    /* Equation 11 */
    deltaHPrime = 2.0 * sqrt(CPrimeProduct) * sin(deltahPrime / 2.0);
    
    /*
     * Step 3
     */
    /* Equation 12 */
    barLPrime = (lab1->l + lab2->l) / 2.0;
    /* Equation 13 */
    barCPrime = (CPrime1 + CPrime2) / 2.0;
    /* Equation 14 */
    hPrimeSum = hPrime1 + hPrime2;
    if (CPrime1 * CPrime2 == 0) {
        barhPrime = hPrimeSum;
    } else {
        if (fabs(hPrime1 - hPrime2) <= deg180InRad)
            barhPrime = hPrimeSum / 2.0;
        else {
            if (hPrimeSum < deg360InRad)
                barhPrime = (hPrimeSum + deg360InRad) / 2.0;
            else
                barhPrime = (hPrimeSum - deg360InRad) / 2.0;
        }
    }
    /* Equation 15 */
    T = 1.0 - (0.17 * cos(barhPrime - deg2Rad(30.0))) +
    (0.24 * cos(2.0 * barhPrime)) +
    (0.32 * cos((3.0 * barhPrime) + deg2Rad(6.0))) -
    (0.20 * cos((4.0 * barhPrime) - deg2Rad(63.0)));
    /* Equation 16 */
    deltaTheta = deg2Rad(30.0) *
    exp(-pow((barhPrime - deg2Rad(275.0)) / deg2Rad(25.0), 2.0));
    /* Equation 17 */
    R_C = 2.0 * sqrt(pow(barCPrime, 7.0) /
                     (pow(barCPrime, 7.0) + pow25To7));
    /* Equation 18 */
    S_L = 1 + ((0.015 * pow(barLPrime - 50.0, 2.0)) /
               sqrt(20 + pow(barLPrime - 50.0, 2.0)));
    /* Equation 19 */
    S_C = 1 + (0.045 * barCPrime);
    /* Equation 20 */
    S_H = 1 + (0.015 * barCPrime * T);
    /* Equation 21 */
    R_T = (-sin(2.0 * deltaTheta)) * R_C;
    
    /* Equation 22 */
    deltaE = sqrt(
                  pow(deltaLPrime / (k_L * S_L), 2.0) +
                  pow(deltaCPrime / (k_C * S_C), 2.0) +
                  pow(deltaHPrime / (k_H * S_H), 2.0) +
                  (R_T * (deltaCPrime / (k_C * S_C)) * (deltaHPrime / (k_H * S_H))));
    
    return (deltaE);
}

/* End of external code */

/*
 ** My comparison function
 */
double comparison(unsigned char *a, unsigned char *b)
{
    LAB c;
    LAB d;
    
    rgb2lab(a, &c);
    rgb2lab(b, &d);
    return CIEDE2000(&c, &d);
}

/*
 ** Generate db assembler data
 */
void generate_db(FILE *output, unsigned char *data, int width, int length, int compress)
{
    int c;
    
    if (compress) {
        pletter(data, length, &data, &length);
        width = 8;
    }
    for (c = 0; c < length; c++) {
        if (c % width == 0)
            fprintf(output, "\tdb ");
        else
            fprintf(output, ",");
        fprintf(output, "$%02x", data[c]);
        if (c % width == width - 1 || (c + 1) == length)
            fprintf(output, "\n");
    }
    if (compress) {
        free(data);
    }
}

/*
 ** Converts a byte to a binary string
 */
char *binary(int data)
{
    static char string[9];
    
    if (use_bitmap == 2) {
        string[0] = (data & 0x80) ? 'X' : '.';
        string[1] = (data & 0x40) ? 'X' : '.';
        string[2] = (data & 0x20) ? 'X' : '.';
        string[3] = (data & 0x10) ? 'X' : '.';
        string[4] = (data & 0x08) ? 'X' : '.';
        string[5] = (data & 0x04) ? 'X' : '.';
        string[6] = (data & 0x02) ? 'X' : '.';
        string[7] = (data & 0x01) ? 'X' : '.';
    } else {
        string[0] = (data & 0x80) ? '1' : '0';
        string[1] = (data & 0x40) ? '1' : '0';
        string[2] = (data & 0x20) ? '1' : '0';
        string[3] = (data & 0x10) ? '1' : '0';
        string[4] = (data & 0x08) ? '1' : '0';
        string[5] = (data & 0x04) ? '1' : '0';
        string[6] = (data & 0x02) ? '1' : '0';
        string[7] = (data & 0x01) ? '1' : '0';
    }
    string[8] = '\0';
    return string;
}

/*
 ** Converts a bitmap to a string
 */
char *nes_color(unsigned char *data)
{
    static char string[9];
    int color[8];
    int c;
    
    for (c = 0; c < 8; c++) {
        string[c] = ((data[0] & (0x80 >> c)) ? 1 : 0)
                  | ((data[8] & (0x80 >> c)) ? 2 : 0);
        string[c] += 0x30;
    }
    string[8] = '\0';
    return string;
}

/*
 ** Converts a bitmap to a string
 */
char *sms_color(unsigned char *data)
{
    static char string[9];
    int color[8];
    int c;
    
    for (c = 0; c < 8; c++) {
        string[c] = ((data[0] & (0x80 >> c)) ? 1 : 0)
                  | ((data[1] & (0x80 >> c)) ? 2 : 0)
                  | ((data[2] & (0x80 >> c)) ? 4 : 0)
                  | ((data[3] & (0x80 >> c)) ? 8 : 0);
        if (string[c] > 9)
            string[c] += 7;
        string[c] += 0x30;
    }
    string[8] = '\0';
    return string;
}

/*
 ** Generate DATA data
 */
void generate_data(FILE *output, unsigned char *data, int width, int length, int compress, int literal)
{
    int c;
    
    if (compress) {
        width = 8;
        pletter(data, length, &data, &length);
        for (c = 0; c < length; c++) {
            if (c % width == 0)
                fprintf(output, "\tDATA BYTE ");
            else
                fprintf(output, ",");
            fprintf(output, "$%02x", data[c]);
            if (c % width == width - 1 || (c + 1) == length)
                fprintf(output, "\n");
        }
        free(data);
    } else if (use_bitmap && !literal) {
        if (nes_mode) {
            for (c = 0; c < length; c += 16) {
                fprintf(output, "\tBITMAP \"%s\"\n", nes_color(&data[c + 0]));
                fprintf(output, "\tBITMAP \"%s\"\n", nes_color(&data[c + 1]));
                fprintf(output, "\tBITMAP \"%s\"\n", nes_color(&data[c + 2]));
                fprintf(output, "\tBITMAP \"%s\"\n", nes_color(&data[c + 3]));
                fprintf(output, "\tBITMAP \"%s\"\n", nes_color(&data[c + 4]));
                fprintf(output, "\tBITMAP \"%s\"\n", nes_color(&data[c + 5]));
                fprintf(output, "\tBITMAP \"%s\"\n", nes_color(&data[c + 6]));
                fprintf(output, "\tBITMAP \"%s\"\n\n", nes_color(&data[c + 7]));
            }
        } else if (sms_mode) {
            for (c = 0; c < length; c += 32) {
                fprintf(output, "\tBITMAP \"%s\"\n", sms_color(&data[c + 0]));
                fprintf(output, "\tBITMAP \"%s\"\n", sms_color(&data[c + 4]));
                fprintf(output, "\tBITMAP \"%s\"\n", sms_color(&data[c + 8]));
                fprintf(output, "\tBITMAP \"%s\"\n", sms_color(&data[c + 12]));
                fprintf(output, "\tBITMAP \"%s\"\n", sms_color(&data[c + 16]));
                fprintf(output, "\tBITMAP \"%s\"\n", sms_color(&data[c + 20]));
                fprintf(output, "\tBITMAP \"%s\"\n", sms_color(&data[c + 24]));
                fprintf(output, "\tBITMAP \"%s\"\n\n", sms_color(&data[c + 28]));
            }
        } else {
            for (c = 0; c < length; c += 8) {
                fprintf(output, "\tBITMAP \"%s\"\n", binary(data[c + 0]));
                fprintf(output, "\tBITMAP \"%s\"\n", binary(data[c + 1]));
                fprintf(output, "\tBITMAP \"%s\"\n", binary(data[c + 2]));
                fprintf(output, "\tBITMAP \"%s\"\n", binary(data[c + 3]));
                fprintf(output, "\tBITMAP \"%s\"\n", binary(data[c + 4]));
                fprintf(output, "\tBITMAP \"%s\"\n", binary(data[c + 5]));
                fprintf(output, "\tBITMAP \"%s\"\n", binary(data[c + 6]));
                fprintf(output, "\tBITMAP \"%s\"\n\n", binary(data[c + 7]));
            }
        }
    } else {
        for (c = 0; c < length; c++) {
            if (c % width == 0)
                fprintf(output, "\tDATA BYTE ");
            else
                fprintf(output, ",");
            fprintf(output, "$%02x", data[c]);
            if (c % width == width - 1 || (c + 1) == length)
                fprintf(output, "\n");
        }
    }
}

/*
 ** Main program
 */
int main(int argc, char *argv[])
{
    FILE *a;
    FILE *output;
    char *p;
    int x;
    int y;
    int c;
    int d;
    int e;
    int n;
    unsigned char buffer[54 + 1024];    /* Header and palette */
    unsigned char color_replacement[32];
    int color_replacement_size = 0;
    unsigned char mapping[16];
    int b;
    int g;
    int r;
    int offset;
    int max1;
    int color1;
    int max2;
    int color2;
    int x1;
    int y1;
    int x2;
    int sig_sprite;
    int triple_color;
    int inline_sprites[192];
    
    int png_file;
    int bmp_format;
    
    int arg;
    int bad;
    int sprite_mode = 0;
    int flip_x = 0;
    int flip_y = 0;
    int magic_sprites = 0;
    int photo = 0;
    int tiled = 0;
    int start_tile = 0;
    int cvbasic = 0;
    int remove_stub = 0;
    int pletter = 0;
    int direct = 0;
    char *output_file = NULL;
    char *label;
    
    time_t actual;
    struct tm *date;
    
    actual = time(0);
    date = localtime(&actual);
    if (argc < 3) {
        fprintf(stderr, "\nTMSColor: Converter from BMP/PNG to TMS9928 format\n");
        fprintf(stderr, VERSION "  by Oscar Toledo G. http://nanochess.org\n");
        fprintf(stderr, "\n");
        fprintf(stderr, "Usage:\n\n");
        fprintf(stderr, "    tmscolor [options] image.bmp image.asm [label]\n");
        fprintf(stderr, "        Creates image for use with assembler code\n\n");
        fprintf(stderr, "    -sms   Sega Master System mode.\n");
        fprintf(stderr, "           Allows up to 448 tiles using -t\n");
        fprintf(stderr, "    -nes   NES/Famicom mode.\n");
        fprintf(stderr, "           Allows up to 256 tiles using -t\n");
        fprintf(stderr, "    -b     Generates CVBasic source code.\n");
        fprintf(stderr, "    -n     Removes CVBasic stub code for displaying.\n");
        fprintf(stderr, "    -z     Output file is compressed with Pletter.\n");
        fprintf(stderr, "    -s     Process tiles in chunks of 16 pixels high (sprites).\n");
        fprintf(stderr, "    -sb    Same as above but generates BITMAP statements.\n");
        fprintf(stderr, "    -t     Generates minimum of tiles required.\n");
        fprintf(stderr, "    -t1    Same but starting at tile 1 (0-255).\n");
        fprintf(stderr, "    -e45d2 Replaces color 4 with 5 and d with 2 before processing.\n");
        fprintf(stderr, "    -fx    Flip image along the X-coordinate (mirror)\n");
        fprintf(stderr, "    -fy    Flip image along the Y-coordinate\n");
        fprintf(stderr, "    -i     Generates BITMAP statements instead of DATA\n");
        fprintf(stderr, "    -i2    Generates BITMAP statements using X and .\n");
        fprintf(stderr, "    -m     Generates magic sprites for areas with more than 2 colors\n");
        fprintf(stderr, "    -p1    Searches best color combination for photo (slow)\n");
        fprintf(stderr, "    -p2    Searches best color combination for photo (2x2 dither) (slow)\n");
        fprintf(stderr, "    -o result.bmp\n");
        fprintf(stderr, "           Outputs the final image, plus highlight of errors (if any).\n");
        fprintf(stderr, "    -d     Direct copy of binary input file to output (can compress).\n");
        fprintf(stderr, "           Useful for getting binary data in ASM or CVBasic code.\n");
        fprintf(stderr, "\n");
        fprintf(stderr, "    Best photo conversion is generated by this command line:\n");
        fprintf(stderr, "      tmscolor -p2 photo.bmp photo.bin\n");
        fprintf(stderr, "    Photos will look better if the contrast is good.\n");
        fprintf(stderr, "    Magic sprites will work only with an image of 256x192 pixels.\n");
        fprintf(stderr, "\n");
        fprintf(stderr, "Thanks to LeandroCorreia for ideas of color conversion.\n");
        fprintf(stderr, "\n");
        return 0;
    }
    arg = 1;
    while (arg < argc && argv[arg][0] == '-') {
        bad = 0;
        c = tolower(argv[arg][1]);
        if (c == 'b') {
            cvbasic = 1;
        } else if (c == 'n') {
            if (tolower(argv[arg][2]) == 'e' && tolower(argv[arg][3]) == 's') {
                nes_mode = 1;
            } else {
                remove_stub = 1;
            }
        } else if (c == 'z') {
            pletter = 1;
        } else if (c == 'f') {
            d = tolower(argv[arg][2]);
            if (d == 'x')
                flip_x = 1;
            else if (d == 'y')
                flip_y = 1;
            else
                bad = 1;
        } else if (c == 's') {     /* -s Sprite mode */
            if (tolower(argv[arg][2]) == 'm' && tolower(argv[arg][3]) == 's') {
                sms_mode = 1;
            } else {
                sprite_mode = 1;
                if (tolower(argv[arg][2]) == 'b')
                    sprite_mode = 2;
            }
        } else if (c == 'e') {  /* -e Color replacement */
            char *ap1 = &argv[arg][2];
            
            while (isxdigit(ap1[0]) && isxdigit(ap1[1])) {
                if (color_replacement_size < 32) {
                    color_replacement[color_replacement_size++] = from_hex(ap1[0]);
                    color_replacement[color_replacement_size++] = from_hex(ap1[1]);
                } else {
                    fprintf(stderr, "Error: too many color replacements in option -e\n");
                }
                ap1 += 2;
            }
        } else if (c == 'm') {  /* -m Magic sprites */
            magic_sprites = 1;
        } else if (c == 'p') {  /* -p Photo */
            d = tolower(argv[arg][2]);
            if (d == '1')
                photo = 1;
            else if (d == '2')
                photo = 2;
            else
                bad = 1;
        } else if (c == 't') {  /* -t Tiled, -t1 Tiles start at 1 */
            char *ap1 = &argv[arg][2];
            
            if (isdigit(ap1[0])) {
                start_tile = atoi(ap1);
                tiled = 1;
            } else if (ap1[0] == '\0') {
                tiled = 1;
            } else {
                bad = 1;
            }
        } else if (c == 'o') {  /* -o file.bmp */
            if (arg + 1 >= argc)
                bad = 1;
            else
                output_file = argv[++arg];
        } else if (c == 'd') {  /* -d */
            direct = 1;
        } else if (c == 'i') {  /* -i Use BITMAP instead of DATA (CVBasic) */
            if (argv[arg][2])
                use_bitmap = atoi(argv[arg] + 2);
            else
                use_bitmap = 1;
            if (use_bitmap == 0)
                use_bitmap = 1;
        } else {
            bad = 1;
        }
        if (bad)
            fprintf(stderr, "Error: Bad option %s\n", argv[arg]);
        arg++;
    }
    if (arg >= argc) {
        fprintf(stderr, "Error: Missing input file name\n");
        exit(2);
    }
    if (nes_mode) {
        if (magic_sprites) {
            fprintf(stderr, "Warning: NES/Famicom mode doesn't require magic sprites\n");
            magic_sprites = 0;
        }
        if (pletter) {
            fprintf(stderr, "Warning: NES/Famicom doesn't support Pletter compression\n");
            pletter = 0;
        }
        if (cvbasic == 0) {
            fprintf(stderr, "Warning: NES/Famicom only can generate CVBasic code\n");
            cvbasic = 1;
        }
        use_bitmap = 1;
        if (start_tile < 0) {
            fprintf(stderr, "Warning: Negative start tile, using 0\n");
            start_tile = 0;
        }
        if (start_tile > 255) {
            fprintf(stderr, "Warning: Out of range start tile, using 255\n");
            start_tile = 255;
        }
    } else if (sms_mode) {
        if (magic_sprites) {
            fprintf(stderr, "Warning: Sega Master System mode doesn't require magic sprites\n");
            magic_sprites = 0;
        }
        if (start_tile < 0) {
            fprintf(stderr, "Warning: Negative start tile, using 0\n");
            start_tile = 0;
        }
        if (start_tile > 447) {
            fprintf(stderr, "Warning: Out of range start tile, using 447\n");
            start_tile = 447;
        }
    } else {
        if (start_tile < 0) {
            fprintf(stderr, "Warning: Negative start tile, using 0\n");
            start_tile = 0;
        }
        if (start_tile > 255) {
            fprintf(stderr, "Warning: Out of range start tile, using 255\n");
            start_tile = 255;
        }
    }
    /* Generate flip bit table */
    for (c = 0; c < 256; c++) {
        flip_bits[c] = ((c & 0x80) >> 7) | ((c & 0x40) >> 5) | ((c & 0x20) >> 3) | ((c & 0x10) >> 1) | ((c & 0x08) << 1) | ((c & 0x04) >> 3) | ((c & 0x02) >> 5) | ((c & 0x01) << 7);
    }
    fprintf(stdout, "Processing: %s\n", argv[arg]);
    png_file = 0;
    if (!direct && strlen(argv[arg]) > 4) {
        p = argv[arg] + strlen(argv[arg]);
        if (p[-4] == '.' && tolower(p[-3]) == 'p' && tolower(p[-2]) == 'n' && tolower(p[-1]) == 'g') {
            png_file = 1;
        }
    }
    if (png_file == 0) {
        a = fopen(argv[arg], "rb");
        arg++;
        if (a == NULL) {
            fprintf(stderr, "Error: Missing input file\n");
            exit(2);
        }
        
        /*
         ** Direct filter (for converting binary files or compressing binary files)
         */
        if (direct) {
            unsigned char *p;
            
            fseek(a, 0, SEEK_END);
            c = ftell(a);
            fseek(a, 0, SEEK_SET);
            bitmap = malloc(c);
            if (bitmap == NULL) {
                fprintf(stderr, "Error: Not enough memory to read file\n");
                exit(2);
            }
            fread(bitmap, 1, c, a);
            p = bitmap + c;
            fclose(a);
            if (arg >= argc) {
                fprintf(stderr, "Error: Missing output file name\n");
                free(bitmap);
                exit(2);
            }
            output = fopen(argv[arg], "w");
            if (output == NULL) {
                fprintf(stderr, "Error: Couldn't open output file '%s'\n", argv[arg]);
                free(bitmap);
                exit(2);
            }
            arg++;
            if (arg < argc) {
                label = argv[arg];
            } else {
                label = "binary";
            }
            if (cvbasic) {
                fprintf(output, "\t' TMSColor " VERSION "\n");
                fprintf(output, "\t' Command: ");
                for (c = 0; c < argc; c++) {
                    char *b;
                    
                    b = strchr(argv[c], ' ');
                    if (b != NULL)
                        fprintf(output, "\"%s\" ", argv[c]);
                    else
                        fprintf(output, "%s ", argv[c]);
                }
                fprintf(output, "\n");
                fprintf(output, "\t' Created: %s\n", asctime(date));
            } else {
                fprintf(output, "\t; TMSColor " VERSION "\n");
                fprintf(output, "\t; Command: ");
                for (c = 0; c < argc; c++) {
                    char *b;
                    
                    b = strchr(argv[c], ' ');
                    if (b != NULL)
                        fprintf(output, "\"%s\" ", argv[c]);
                    else
                        fprintf(output, "%s ", argv[c]);
                }
                fprintf(output, "\n");
                fprintf(output, "\t; Created: %s\n", asctime(date));
            }
            fprintf(output, "%s:\n", label);
            if (cvbasic)
                generate_data(output, bitmap, 8, p - bitmap, pletter, 1);
            else
                generate_db(output, bitmap, 8, p - bitmap, pletter);
            fclose(output);
            exit(0);
        }
        
        /*
         ** Read the BMP file
         */
        fread(buffer, 1, 54 + 1024, a);
        if (buffer[0] != 'B' || buffer[1] != 'M') {
            fprintf(stderr, "Error: The input file is not in BMP format\n");
            fclose(a);
            exit(3);
        }
        bmp_format = buffer[0x1c];
        if (bmp_format != 8 && bmp_format != 24 && bmp_format != 32) {
            fprintf(stderr, "Error: The input file is in %d bits format (not 8/24/32)\n", bmp_format);
            fclose(a);
            exit(3);
        }
        if (bmp_format == 8) {
            if (buffer[0x2e] != 0 || (buffer[0x2f] != 0 && buffer[0x2f] != 1)) {
                fprintf(stderr, "Error: Unsupported palette for 8 bits format\n");
                fclose(a);
                exit(3);
            }
        }
        if (buffer[0x1e] != 0 || buffer[0x1f] != 0 || buffer[0x20] != 0 || buffer[0x21] != 0) {
            fprintf(stderr, "Error: Cannot handle compressed input files (codec 0x%08x)\n", (buffer[0x21] << 24) | (buffer[0x20] << 16) | (buffer[0x1f] << 8) | (buffer[0x1e]));
            fclose(a);
            exit(3);
        }
        size_x = buffer[0x12] | (buffer[0x13] << 8);
        size_y = buffer[0x16] | (buffer[0x17] << 8);
        if (size_y >= 32768)
            size_y -= 65536;
        if (size_y >= 0) {
            n = 0;
        } else {
            size_y = -size_y;
            n = 1;
        }
        if (size_x == 0 || size_y == 0) {
            fprintf(stderr, "Error: There's a weird BMP file in the input. I'm scared...\n");
            fclose(a);
            exit(3);
        }
        image = malloc(size_x * size_y * 3);
        if (image == NULL) {
            fprintf(stderr, "Error: Couldn't allocate memory for image\n");
            fclose(a);
            exit(3);
        }
        /* Read image and normalize as 24-bit RGB */
        fseek(a, buffer[10] | (buffer[11] << 8) | (buffer[12] << 16) | (buffer[13] << 24), SEEK_SET);
        for (y = n ? 0 : size_y - 1; n ? y < size_y : y >= 0; n ? y++ : y--) {
            for (x = 0; x < size_x; x++) {
                if (bmp_format == 8) {            /* 256 color */
                    fread(buffer, 1, 1, a);
                    memcpy(buffer, buffer + 54 + buffer[0] * 4, 4);
                } else if (bmp_format == 24) {    /* 24 bits */
                    fread(buffer, 1, 3, a);
                } else {                            /* 32 bits */
                    fread(buffer, 1, 4, a);
                }
                image[(y * size_x + x) * 3 + 0] = buffer[2];
                image[(y * size_x + x) * 3 + 1] = buffer[1];
                image[(y * size_x + x) * 3 + 2] = buffer[0];
            }
        }
        fclose(a);
    } else {
        image = NULL;
        size_x = 0;
        size_y = 0;
        c = lodepng_decode24_file(&image, (unsigned *) &size_x, (unsigned *) &size_y, argv[arg]);
        if (c) {
            fprintf(stderr, "Error: %s\n", lodepng_error_text(c));
            exit(2);
        }
        arg++;
    }
    if ((size_x & 7) != 0) {
        fprintf(stderr, "Error: The input file doesn't measure a multiple of 8 in X size (it's %d pixels)\n", size_x);
        exit(3);
    }
    if ((size_y & 7) != 0) {
        fprintf(stderr, "Error: The input file doesn't measure a multiple of 8 in Y size (it's %d pixels)\n", size_y);
        exit(3);
    }

    if (nes_mode) {
        bitmap = malloc(size_x * size_y / 4);
        color = NULL;
        pattern = malloc((size_x + 7) / 8 * size_y / 8);
    } else if (sms_mode) {
        bitmap = malloc(size_x * size_y / 2);
        color = NULL;
        pattern = malloc((size_x + 7) / 8 * size_y / 8 * 2);    /* Tile is word */
    } else {
        bitmap = malloc(size_x * size_y / 8);
        color = malloc(size_x * size_y / 8);
        pattern = malloc((size_x + 7) / 8 * size_y / 8);
        if (color == NULL) {
            fprintf(stderr, "Error: Unable to allocate memory for color data\n");
            fclose(a);
            exit(3);
        }
    }
    source = malloc(size_x * size_y);
    source2 = malloc(size_x * size_y);
    ignore = malloc(size_x * size_y);
    if (bitmap == NULL || pattern == NULL || source == NULL || source2 == NULL || ignore == NULL) {
        fprintf(stderr, "Error: Unable to allocate memory for bitmap\n");
        fclose(a);
        exit(3);
    }
    for (y = 0; y < size_y; y++) {
        
        /*
         ** If asked for photo quality...
         */
        if (photo && nes_mode) {
            double g;
            
            /*
             ** Normal image
             */
            for (x = 0; x < size_x; x++) {
                int best_color;
                
                buffer[0] = image[(y * size_x + x) * 3 + 2];
                buffer[1] = image[(y * size_x + x) * 3 + 1];
                buffer[2] = image[(y * size_x + x) * 3 + 0];
                best_color = ((buffer[0] * 30 + buffer[1] * 59 + buffer[2] * 11) / 100) >> 6;
                source[(flip_y ? size_y - 1 - y : y) * size_x + (flip_x ? size_x - 1 - x : x)] = best_color;
            }
        } else if (photo && sms_mode) {
            double g;
            
            /*
             ** Normal image
             */
            for (x = 0; x < size_x; x++) {
                int best_color;
                double best_difference;
                
                buffer[0] = image[(y * size_x + x) * 3 + 2];
                buffer[1] = image[(y * size_x + x) * 3 + 1];
                buffer[2] = image[(y * size_x + x) * 3 + 0];
                best_color = 0;
                best_difference = 1e38;
                for (c = 0; c < 16; c++) {
                    g = comparison(&buffer[0], &colors[c * 3 + 0]);
                    if (g < best_difference) {
                        best_difference = g;
                        best_color = c;
                    }
                }
                source[(flip_y ? size_y - 1 - y : y) * size_x + (flip_x ? size_x - 1 - x : x)] = best_color;
            }
        } else if (photo) {
            int best_combo;
            double best_combo_difference;
            int best_combo_dither;
            double combo_difference;
            int combo_dither;
            double g;
            double r;
            double b;
            
            for (x = 0; x < size_x; ) {
                for (e = 0; e < 32; e += 4) {
                    buffer[e + 0] = image[(y * size_x + x) * 3 + 2];
                    buffer[e + 1] = image[(y * size_x + x) * 3 + 1];
                    buffer[e + 2] = image[(y * size_x + x) * 3 + 0];
                    x++;
                }
                x -= 8;
                best_combo_difference = 1e38;
                
                /*
                 ** Try all possible combinations of two colors (15 * 15 = 225 possibilities)
                 */
                for (c = 1; c < 16; c++) {
                    for (d = 1; d < 16; d++) {
                        if (d <= c)
                            continue;
                        combo_difference = 0;
                        combo_dither = 0;
                        for (e = 0; e < 32; e += 4) {
                            g = comparison(&buffer[e + 0], &colors[c * 3 + 0]);
                            r = comparison(&buffer[e + 0], &colors[d * 3 + 0]);
                            if (photo == 2) {
                                g = (g < r) ? g : r;
                                r = sqrt(sqr(colors[d * 3 + 2] - colors[c * 3 + 2])
                                         + sqr(colors[d * 3 + 1] - colors[d * 3 + 1])
                                         + sqr(colors[d * 3 + 0] - colors[d * 3 + 0]));
                                if (r < 96) {
                                    unsigned char extra[3];
                                    
                                    extra[0] = (colors[d * 3 + 0] + colors[c * 3 + 0]) / 2;
                                    extra[1] = (colors[d * 3 + 1] + colors[c * 3 + 1]) / 2;
                                    extra[2] = (colors[d * 3 + 2] + colors[c * 3 + 2]) / 2;
                                    
                                    r = comparison(&buffer[e + 0], &extra[0]);
                                    g = (g < r) ? g : r;
                                    combo_dither = 1;
                                }
                                combo_difference += g;
                            } else {
                                combo_difference += (g < r) ? g : r;
                            }
                        }
                        if (combo_difference < best_combo_difference) {
                            best_combo = (c << 4) | d;
                            best_combo_difference = combo_difference;
                            best_combo_dither = combo_dither;
                        }
                    }
                }
                for (c = 0; c < color_replacement_size; c += 2) {
                    if ((best_combo >> 4) == color_replacement[c]) {
                        best_combo = (best_combo & 0x0f) | (color_replacement[c + 1] << 4);
                        break;
                    }
                }
                for (c = 0; c < color_replacement_size; c += 2) {
                    if ((best_combo & 0x0f) == color_replacement[c]) {
                        best_combo = (best_combo & 0xf0) | color_replacement[c + 1];
                        break;
                    }
                }
                c = (best_combo >> 4) & 0x0f;
                d = best_combo & 0x0f;
                for (e = 0; e < 32; e += 4) {
                    g = comparison(&buffer[e + 0], &colors[c * 3 + 0]);
                    r = comparison(&buffer[e + 0], &colors[d * 3 + 0]);
                    if (best_combo_dither) {
                        unsigned char extra[3];
                        
                        extra[0] = (colors[d * 3 + 0] + colors[c * 3 + 0]) / 2;
                        extra[1] = (colors[d * 3 + 1] + colors[c * 3 + 1]) / 2;
                        extra[2] = (colors[d * 3 + 2] + colors[c * 3 + 2]) / 2;
                        
                        b = comparison(&buffer[e + 0], &extra[0]);
                        if (b < g && b < r) {
                            source[(flip_y ? size_y - 1 - y : y) * size_x + (flip_x ? size_x - 1 - x : x)] = (y ^ x) & 1 ? c : d;
                        } else {
                            source[(flip_y ? size_y - 1 - y : y) * size_x + (flip_x ? size_x - 1 - x : x)] = (g < r) ? c : d;
                        }
                    } else {
                        source[(flip_y ? size_y - 1 - y : y) * size_x + (flip_x ? size_x - 1 - x : x)] = (g < r) ? c : d;
                    }
                    x++;
                }
            }
        } else {
            
            /*
             ** Normal image
             */
            for (x = 0; x < size_x; x++) {
                int best_color;
                int best_difference;
                
                buffer[0] = image[(y * size_x + x) * 3 + 2];
                buffer[1] = image[(y * size_x + x) * 3 + 1];
                buffer[2] = image[(y * size_x + x) * 3 + 0];
                best_color = 0;
                best_difference = 100000;
                for (c = 0; c < 16; c++) {
                    d = (buffer[2] - colors[c * 3 + 2]) * (buffer[2] - colors[c * 3 + 2])
                    + (buffer[1] - colors[c * 3 + 1]) * (buffer[1] - colors[c * 3 + 1])
                    + (buffer[0] - colors[c * 3 + 0]) * (buffer[0] - colors[c * 3 + 0]);
                    if (d < best_difference) {
                        best_difference = d;
                        best_color = c;
                    }
                }
                for (c = 0; c < color_replacement_size; c += 2) {
                    if (best_color == color_replacement[c]) {
                        best_color = color_replacement[c + 1];
                        break;
                    }
                }
                source[(flip_y ? size_y - 1 - y : y) * size_x + (flip_x ? size_x - 1 - x : x)] = best_color;
            }
        }
    }

    if (size_x == 256 && size_y == 192 && magic_sprites == 1) {
        memset(inline_sprites, 0, sizeof(inline_sprites));
        memset(sprites, 0, sizeof(sprites));
        memset(attr, 0xd1, sizeof(attr));
        sig_sprite = 0;
        
        do {
        hack1:
            triple_color = check_triple_color();
            
            
            if (triple_color == 0)
                break;
            
            //
            // Do passes over image replacing more than 2 colors with
            // 16x16 sprites
            //
            for (c = 0; c < 192; c++) {
                for (d = 0; d < 256; d += 8) {
                    if (ignore[c * size_x + d])
                        continue;
                    //                fprintf(stderr, "Generating sprite at %d,%d (%d,%d,%d)\n", d * 8, c * 8, usage[c][d][0], usage[c][d][1], usage[c][d][2]);
                    
                    b = 0;
                    color2 = 0;
                    max2 = -1;
                    while (1) {
                        g = usage[c / 8][d / 8][b];
                        if (g == 255)
                            break;
                        y1 = c;
                        x1 = d;
                        for (x = 0; x < 16; x++) {
                            for (y = 0; y < 16; y++) {
                                if (x + x1 < 256 && y + y1 < 192 && source[(y + y1) * size_x + x + x1] == g && !ignore[(y + y1) * size_x + x + x1])
                                    break;
                            }
                            if (y < 16) {
                                x1 += x;
                                break;
                            }
                        }
                        x2 = 16;
                        for (x = 15; x >= 0; x--) {
                            for (y = 0; y < 16; y++) {
                                if (x + x1 < 256 && y + y1 < 192 && source[(y + y1) * size_x + x + x1] == g && !ignore[(y + y1) * size_x + x + x1])
                                    break;
                            }
                            if (y < 16) {
                                x2 = x + 1;
                                break;
                            }
                        }
                        for (x = x2 - 16; x < 0; x++) {
                            for (y = 0; y < 16; y++) {
                                if (x + x1 < 256 && y + y1 < 192 && source[(y + y1) * size_x + x + x1] == g && !ignore[(y + y1) * size_x + x + x1])
                                    break;
                            }
                            if (y < 16) {
                                x1 += x;
                                break;
                            }
                        }
                        for (y = 15; y >= 0; y--) {
                            for (x = 0; x < 16; x++) {
                                if (x + x1 < 256 && y + y1 < 192 && source[(y + y1) * size_x + x + x1] == g && !ignore[(y + y1) * size_x + x + x1])
                                    break;
                            }
                            if (x < 16) {
                                y1 -= 15 - y;
                                if (y1 < 0)
                                    y1 = 0;
                                break;
                            }
                        }
                        memcpy(source2, source,size_x * size_y);
                        for (y = 0; y < 16; y++) {
                            for (x = 0; x < 16; x++) {
                                if (x == 0 || ((x + x1) & 7) == 0 && x + x1 < 256 && y + y1 < 192) {
                                    color1 = 255;
                                    e = (x1 + x) & 0xf8;
                                    do {
                                        if (source[(y + y1) * size_x + e] != g)
                                            color1 = source[(y + y1) * size_x + e];
                                    } while (++e & 7) ;
                                    if (color1 == 255)
                                        color1 = 1;
                                }
                                if (y + y1 < 192 && x + x1 < 256) {
                                    if (source[(y + y1) * size_x + x + x1] == g && !ignore[(y + y1) * size_x + x + x1]) {
                                        source[(y + y1) * size_x + x + x1] = color1;
                                    }
                                }
                            }
                        }
                        max1 = triple_color - check_triple_color();
                        memcpy(source, source2, size_x * size_y);
                        check_triple_color();
                        if (max1 > max2) {
                            max2 = max1;
                            color2 = g;
                        }
                        b++;
                    }
                    if (sig_sprite >= 32) {
                        fprintf(stderr, "Error: More than 32 sprites\n");
                        goto hack;
                    }
                    y1 = c;
                    x1 = d;
                    for (x = 0; x < 16; x++) {
                        for (y = 0; y < 16; y++) {
                            if (x + x1 < 256 && y + y1 < 192 && source[(y + y1) * size_x + x + x1] == color2 && !ignore[(y + y1) * size_x + x + x1])
                                break;
                        }
                        if (y < 16) {
                            x1 += x;
                            break;
                        }
                    }
                    x2 = 16;
                    for (x = 15; x >= 0; x--) {
                        for (y = 0; y < 16; y++) {
                            if (x + x1 < 256 && y + y1 < 192 && source[(y + y1) * size_x + x + x1] == color2 && !ignore[(y + y1) * size_x + x + x1])
                                break;
                        }
                        if (y < 16) {
                            x2 = x + 1;
                            break;
                        }
                    }
                    for (x = x2 - 16; x < 0; x++) {
                        for (y = 0; y < 16; y++) {
                            if (x + x1 < 256 && y + y1 < 192 && source[(y + y1) * size_x + x + x1] == color2 && !ignore[(y + y1) * size_x + x + x1])
                                break;
                        }
                        if (y < 16) {
                            x1 += x;
                            break;
                        }
                    }
                    for (y = 15; y >= 0; y--) {
                        for (x = 0; x < 16; x++) {
                            if (x + x1 < 256 && y + y1 < 192 && source[(y + y1) * size_x + x + x1] == color2 && !ignore[(y + y1) * size_x + x + x1])
                                break;
                        }
                        if (x < 16) {
                            y1 -= 15 - y;
                            if (y1 < 0)
                                y1 = 0;
                            break;
                        }
                    }
                    g = sig_sprite * 4;
                    attr[g] = y1 - 1;
                    attr[g + 1] = x1;
                    attr[g + 2] = sig_sprite * 4;
                    attr[g + 3] = color2;
                    //                fprintf(stderr, "%02x,%02x,%02x,%02x\n", y1 - 1, x1, sig_sprite * 4, color2);
                    g = sig_sprite * 32;
                    for (y = 0; y < 16; y++) {
                        for (x = 0; x < 16; x++) {
                            if (x == 0 || ((x + x1) & 7) == 0 && x + x1 < 256 && y + y1 < 192) {
                                color1 = 255;
                                e = (x1 + x) & 0xf8;
                                do {
                                    if (source[(y + y1) * size_x + e] != color2)
                                        color1 = source[(y + y1) * size_x + e];
                                } while (++e & 7) ;
                                if (color1 == 255)
                                    color1 = 1;
                            }
                            if (y + y1 < 192 && x + x1 < 256) {
                                if (source[(y + y1) * size_x + x + x1] == color2 && !ignore[(y + y1) * size_x + x + x1]) {
                                    sprites[g + y + (x / 8) * 16] |= 0x80 >> (x & 7);
                                    source[(y + y1) * size_x + x + x1] = color1;
                                }
                            }
                        }
                    }
                    sig_sprite++;
                    goto hack1;
                }
            }
        } while (1) ;
    }
hack:
    //    fprintf(stderr, "Processing image...\n");
    if (nes_mode) { /* NES/Famicom mode */
        for (c = 0; c < size_y; c++) {
            for (d = 0; d < size_x; d += 8) {
                unsigned int word;
                
                word = 0;
                for (e = 0; e < 8; e++) {
                    n = source[c * size_x + d + e];
                    if (n & 1)
                        word |= 0x80 >> e;
                    if (n & 2)
                        word |= 0x8000 >> e;
                }
                offset = c / 8 * (size_x * 2) + (c & 7) + (d * 2);
                bitmap[offset] = word;
                bitmap[offset + 8] = word >> 8;
                offset++;
            }
        }
    } else if (sms_mode) { /* Sega Master System mode */
        for (c = 0; c < size_y; c++) {
            offset = c / 8 * (size_x * 4) + (c & 7) * 4;
            for (d = 0; d < size_x; d += 8) {
                unsigned int word;
                
                word = 0;
                for (e = 0; e < 8; e++) {
                    n = source[c * size_x + d + e];
                    if (n & 1)
                        word |= 0x80 >> e;
                    if (n & 2)
                        word |= 0x8000 >> e;
                    if (n & 4)
                        word |= 0x800000 >> e;
                    if (n & 8)
                        word |= 0x80000000 >> e;
                }
                bitmap[offset] = word;
                bitmap[offset + 1] = word >> 8;
                bitmap[offset + 2] = word >> 16;
                bitmap[offset + 3] = word >> 24;
                offset += 32;
            }
        }
    } else {    /* TMS9918 mode */
        for (c = 0; c < size_y; c++) {
            for (d = 0; d < size_x; d += 8) {
                offset = c / 8 * size_x + (c & 7) + d;
                color1 = -1;
                color2 = -1;
                for (e = 0; e < 8; e++) {
                    if (source[c * size_x + d + e] == color1) {
                    } else if (source[c * size_x + d + e] == color2) {
                    } else if (color1 == -1) {
                        color1 = source[c * size_x + d + e];
                    } else if (color2 == -1) {
                        color2 = source[c * size_x + d + e];
                    } else {
                        fprintf(stderr, "Error: More than 2 colors in stripe %d,%d (found %d with %d and %d already)\n", d, c, source[c * size_x + d + e], color1, color2);
                        for (e = 0; e < 8; e++)
                            source[c * size_x + d + e] |= 0x10;
                        break;
                    }
                }
                if (color1 == -1)
                    color1 = 1;
                if (color2 == -1) {
                    if (color1 == 1)
                        color2 = 15;
                    else
                        color2 = 1;
                }
                if (color1 < color2) {
                    e = color1;
                    color1 = color2;
                    color2 = e;
                }
                for (e = 0; e < 16; e++)
                    mapping[e] = 0;
                mapping[color1] = 0x80;
                mapping[color2] = 0;
                color[offset] = color1 << 4 | color2;
                r = 0;
                for (e = 0; e < 8; e++)
                    r |= mapping[source[c * size_x + d + e] & 0x0f] >> e;
                bitmap[offset] = r;
                offset += 8;
            }
        }
    }
    
    /*
     ** Generate bitmap file (if requested) to see how the conversion went.
     */
    if (output_file != NULL) {
        FILE *a;
        
        a = fopen(output_file, "wb");
        if (a == NULL) {
            fprintf(stderr, "Error: Unable to write output file \"%s\"\n", output_file);
        } else {
            char header[54];
            
            memset(header, 0, sizeof(header));
            header[0x00] = 'B';     /* Header */
            header[0x01] = 'M';
            c = size_x * size_y * 3 + 54;
            header[0x02] = c;       /* Complete size of file */
            header[0x03] = c >> 8;
            header[0x04] = c >> 16;
            header[0x05] = c >> 24;
            c = 54;
            header[0x0a] = c;       /* Size of header plus palette */
            c = 40;
            header[0x0e] = c;       /* Size of header */
            header[0x12] = size_x;
            header[0x13] = size_x >> 8;
            header[0x16] = size_y;
            header[0x17] = size_y >> 8;
            header[0x1a] = 0x01;    /* 1 plane */
            header[0x1c] = 0x18;    /* 24 bits */
            c = size_x * size_y * 3;
            header[0x22] = c;       /* Complete size of file */
            header[0x23] = c >> 8;
            header[0x24] = c >> 16;
            header[0x25] = c >> 24;
            c = 0x0ec4;             /* 96 dpi */
            header[0x26] = c;       /* X */
            header[0x27] = c >> 8;
            header[0x2a] = c;       /* Y */
            header[0x2b] = c >> 8;
            fwrite(header, 1, sizeof(header), a);
            
            for (y = size_y - 1; y >= 0; y--) {
                for (x = 0; x < size_x; x++) {
                    header[0x00] = colors[(source[y * size_x + x] & 0x0f) * 3];
                    header[0x01] = colors[(source[y * size_x + x] & 0x0f) * 3 + 1];
                    header[0x02] = colors[(source[y * size_x + x] & 0x0f) * 3 + 2];
                    c = (header[0x00] * 30 + header[0x01] * 59 + header[0x02] * 11) / 100;
                    if (source[y * size_x + x] & 0x10) {  /* Error */
                        header[0x00] = c / 8 + 0x80;
                        header[0x01] = c / 4 + 0x80;
                        header[0x02] = c / 2 + 0x80;
                    }
                    fwrite(header, 1, 3, a);
                }
            }
            fclose(a);
        }
    }
    
    /*
     ** Proceed to open the final output file.
     */
    if (arg >= argc) {
        fprintf(stderr, "Error: Missing output file name\n");
        free(bitmap);
        free(color);
        exit(2);
    }
    output = fopen(argv[arg], "w");
    if (output == NULL) {
        fprintf(stderr, "Error: Couldn't open output file '%s'\n", argv[arg]);
        free(bitmap);
        free(color);
        exit(2);
    }
    arg++;
    if (arg < argc) {
        label = argv[arg];
    } else {
        label = "image";
    }

    /*
     ** Headers
     */
    if (cvbasic) {
        fprintf(output, "\t' TMSColor " VERSION "\n");
        fprintf(output, "\t' Command: ");
        for (c = 0; c < argc; c++) {
            char *b;
            
            b = strchr(argv[c], ' ');
            if (b != NULL)
                fprintf(output, "\"%s\" ", argv[c]);
            else
                fprintf(output, "%s ", argv[c]);
        }
        fprintf(output, "\n");
        fprintf(output, "\t' Created: %s\n", asctime(date));
    } else {
        fprintf(output, "\t; TMSColor " VERSION "\n");
        fprintf(output, "\t; Command: ");
        for (c = 0; c < argc; c++) {
            char *b;
            
            b = strchr(argv[c], ' ');
            if (b != NULL)
                fprintf(output, "\"%s\" ", argv[c]);
            else
                fprintf(output, "%s ", argv[c]);
        }
        fprintf(output, "\n");
        fprintf(output, "\t; Created: %s\n", asctime(date));
    }
    
    /*
     ** Tiled mode
     */
    if (tiled) {
        if (nes_mode) { /* NES/Famicom */
            static unsigned char bit[256 * 16];
            int start = start_tile;
            int current = start;
            int total_tiles;
            
            current = start;
            memset(bit, 0, sizeof(bit));
            for (c = 0; c < size_x / 8 * size_y / 8; c++) {
                int flags;
                
                flags = 0;
                for (e = start; e < (current > 256 ? 256 : current); e++) {
                    if (memcmp(&bit[e * 16], &bitmap[c * 16], 16) == 0)
                        break;
                }
                if (e == current) {
                    if (current == 256) {
                        fprintf(stderr, "Error: More than 256 tiles required at %d,%d\n", c % (size_x / 8), c / (size_x / 8));
                        current++;
                    } else if (current > 256) {
                        current++;
                    } else {
                        memcpy(&bit[e * 16], &bitmap[c * 16], 16);
                        current++;
                    }
                }
                pattern[c] = e;
            }
            total_tiles = current - start;
            fprintf(stderr, "Total used tiles: %d ($%02x-$%02x)\n", total_tiles, start_tile, current - 1);
            if (remove_stub)
                fprintf(output, "\t'\n\t' Recommended code:\n");
            else
                fprintf(output, "\t' Display image.\n");
            fprintf(output, "\t%sSCREEN DISABLE\n", remove_stub ? "' " : "");
            fprintf(output, "\t%sSCREEN %s_pattern,0,0,%d,%d,%d\n", remove_stub ? "' " : "", label, size_x / 8, size_y / 8, size_x / 8);
            fprintf(output, "\t%sWAIT\n", remove_stub ? "' " : "");
            fprintf(output, "\t%sSCREEN ENABLE\n", remove_stub ? "' " : "");
            if (remove_stub)
                fprintf(output, "\t'\n");
            else
                fprintf(output, "\tWHILE 1: WEND\n\n");
            fprintf(output, "\t' Start tile = %d. Total_tiles = %d\n", start_tile, total_tiles);
            fprintf(output, "\tCHRROM 0 PATTERN %d\n", start_tile);
            generate_data(output, bit + start_tile * 16, 8, total_tiles * 16, 0, 0);
            fprintf(output, "\n");
            fprintf(output, "\t' Width = %d, height = %d\n", size_x / 8, size_y / 8);
            fprintf(output, "%s_pattern:\n", label);
            generate_data(output, pattern, size_x / 8, size_x / 8 * size_y / 8, 0, 1);
        } else if (sms_mode) { /* Sega Master System */
            static unsigned char bit[512 * 32];
            int start = start_tile;
            int current = start;
            int total_tiles;
            
            current = start;
            memset(bit, 0, sizeof(bit));
            for (c = 0; c < size_x / 8 * size_y / 8; c++) {
                int flags;
            
                flags = 0;
                for (e = start; e < (current > 448 ? 448 : current); e++) {
                    if (memcmp(&bit[e * 32], &bitmap[c * 32], 32) == 0)
                        break;
                    if (check_flip_x(&bit[e * 32], &bitmap[c * 32])) {
                        flags |= 0x02;
                        break;
                    }
                    if (check_flip_y(&bit[e * 32], &bitmap[c * 32])) {
                        flags |= 0x04;
                        break;
                    }
                    if (check_flip_xy(&bit[e * 32], &bitmap[c * 32])) {
                        flags |= 0x06;
                        break;
                    }
                }
                if (e == current) {
                    if (current == 256) {
                        fprintf(stderr, "Warning: Tile 256 being used at %d,%d avoid if using sprites\n", c % (size_x / 8), c / (size_x / 8));
                    }
                    if (current == 448) {
                        fprintf(stderr, "Error: More than 448 tiles required at %d,%d\n", c % (size_x / 8), c / (size_x / 8));
                        current++;
                    } else if (current > 448) {
                        current++;
                    } else {
                        memcpy(&bit[e * 32], &bitmap[c * 32], 32);
                        current++;
                    }
                }
                pattern[c * 2 + 0] = e;
                pattern[c * 2 + 1] = flags | ((e >> 8) & 1);
            }
            total_tiles = current - start;
            fprintf(stderr, "Total used tiles: %d ($%02x-$%02x)\n", total_tiles, start_tile, current - 1);
            if (cvbasic) {
                if (remove_stub)
                    fprintf(output, "\t'\n\t' Recommended code:\n");
                else
                    fprintf(output, "\t' Display image.\n");
                fprintf(output, "\t%sMODE 4\n", remove_stub ? "' " : "");
                if (!pletter && total_tiles > 255) {
                    fprintf(output, "\t%sDEFINE CHAR %s%d,%d,%s_char_0\n", remove_stub ? "' " : "", pletter ? "PLETTER " : "", start_tile, 255, label);
                    fprintf(output, "\t%sDEFINE CHAR %s%d,%d,%s_char_1\n", remove_stub ? "' " : "", pletter ? "PLETTER " : "", start_tile, total_tiles - 255, label);
                } else {
                    fprintf(output, "\t%sDEFINE CHAR %s%d,%d,%s_char\n", remove_stub ? "' " : "", pletter ? "PLETTER " : "", start_tile, total_tiles, label);
                }
                fprintf(output, "\t%sSCREEN %s_pattern,0,0,%d,%d,%d\n", remove_stub ? "' " : "", label, size_x / 8, size_y / 8, size_x / 8);
                if (remove_stub)
                    fprintf(output, "\t'\n");
                else
                    fprintf(output, "\tWHILE 1: WEND\n\n");
                fprintf(output, "\t' Start tile = %d. Total_tiles = %d\n", start_tile, total_tiles);
                if (!pletter && total_tiles > 255) {
                    fprintf(output, "%s_char_0:\n", label);
                    generate_data(output, bit + start_tile * 32, 8, 255 * 32, pletter, 0);
                    fprintf(output, "%s_char_1:\n", label);
                    generate_data(output, bit + (start_tile + 255) * 32, 8, (total_tiles - 255) * 32, pletter, 0);
                } else {
                    fprintf(output, "%s_char:\n", label);
                    generate_data(output, bit + start_tile * 32, 8, total_tiles * 32, pletter, 0);
                }
                fprintf(output, "\n");
                fprintf(output, "\t' Width = %d, height = %d\n", size_x / 8, size_y / 8);
                fprintf(output, "%s_pattern:\n", label);
                generate_data(output, pattern, size_x / 8 * 2, size_x / 8 * 2 * size_y / 8, 0, 1);
            } else {
                fprintf(output, "\t;\n");
                fprintf(output, "\t; Start tile = %d. Total_tiles = %d\n", start_tile, total_tiles);
                fprintf(output, "\t;\n");
                fprintf(output, "%s_char:\n", label);
                generate_db(output, bit + start_tile * 32, 8, total_tiles * 32, pletter);
                fprintf(output, "\n");
                fprintf(output, "\t; Width = %d, height = %d\n", size_x / 8, size_y / 8);
                fprintf(output, "%s_pattern:\n", label);
                generate_db(output, pattern, size_x / 8 * 2, size_x / 8 * 2 * size_y / 8, 0);
            }
        } else {    /* TMS9918 */
            static unsigned char bit[256 * 8];
            static unsigned char col[256 * 8];
            int start = start_tile;
            int current = start;
            int total_tiles;
            
            memset(bit, 0, sizeof(bit));
            memset(col, 0, sizeof(col));
            for (c = 0; c < size_x / 8 * size_y / 8; c++) {
                for (e = start; e < (current > 256 ? 256 : current); e++) {
                    if (memcmp(&bit[e * 8], &bitmap[c * 8], 8) == 0
                        && memcmp(&col[e * 8], &color[c * 8], 8) == 0)
                        break;
                }
                if (e == current) {
                    if (current == 256) {
                        fprintf(stderr, "Error: Too many tiles at %d,%d\n", c % (size_x / 8), c / (size_x / 8));
                        current++;
                    } else if (current > 256) {
                        current++;
                    } else {
                        memcpy(&bit[e * 8], &bitmap[c * 8], 8);
                        memcpy(&col[e * 8], &color[c * 8], 8);
                        current++;
                    }
                }
                pattern[c] = e;
            }
            total_tiles = current - start;
            fprintf(stderr, "Total used tiles: %d ($%02x-$%02x)\n", total_tiles, start_tile, current - 1);
            if (cvbasic) {
                if (remove_stub)
                    fprintf(output, "\t'\n\t' Recommended code:\n");
                else
                    fprintf(output, "\t' Display image.\n");
                fprintf(output, "\t%sMODE 0\n", remove_stub ? "' ": "");
                fprintf(output, "\t%sDEFINE CHAR %s%d,%d,%s_char\n", remove_stub ? "' ": "", pletter ? "PLETTER " : "", start_tile, total_tiles, label);
                fprintf(output, "\t%sDEFINE COLOR %s%d,%d,%s_color\n", remove_stub ? "' ": "", pletter ? "PLETTER " : "", start_tile, total_tiles, label);
                if (magic_sprites)
                    fprintf(output, "\t%sGOSUB %s_show\n", remove_stub ? "' ": "", label);
                fprintf(output, "\t%sSCREEN %s_pattern,0,0,%d,%d,%d\n", remove_stub ? "' ": "", label, size_x / 8, size_y / 8, size_x / 8);
                if (remove_stub)
                    fprintf(output, "\t'\n");
                else
                    fprintf(output, "\tWHILE 1: WEND\n\n");
                fprintf(output, "\t' Start tile = %d. Total_tiles = %d\n", start_tile, total_tiles);
                fprintf(output, "%s_char:\n", label);
                generate_data(output, bit + start_tile * 8, 8, total_tiles * 8, pletter, 0);
                fprintf(output, "\n");
                fprintf(output, "%s_color:\n", label);
                generate_data(output, col + start_tile * 8, 8, total_tiles * 8, pletter, 1);
                fprintf(output, "\n");
                fprintf(output, "\t' Width = %d, height = %d\n", size_x / 8, size_y / 8);
                fprintf(output, "%s_pattern:\n", label);
                generate_data(output, pattern, size_x / 8, size_x / 8 * size_y / 8, 0, 1);
            } else {
                fprintf(output, "\t;\n");
                fprintf(output, "\t; Start tile = %d. Total_tiles = %d\n", start_tile, total_tiles);
                fprintf(output, "\t;\n");
                fprintf(output, "%s_char:\n", label);
                generate_db(output, bit + start_tile * 8, 8, total_tiles * 8, pletter);
                fprintf(output, "\n");
                fprintf(output, "%s_color:\n", label);
                generate_db(output, col + start_tile * 8, 8, total_tiles * 8, pletter);
                fprintf(output, "\n");
                fprintf(output, "\t; Width = %d, height = %d\n", size_x / 8, size_y / 8);
                fprintf(output, "%s_pattern:\n", label);
                generate_db(output, pattern, size_x / 8, size_x / 8 * size_y / 8, 0);
            }
        }
    } else if (sprite_mode) {   /* Bitmap for sprites */
        if (nes_mode) { /* NES/Famicom */
            unsigned char *final_bitmap;
            unsigned char *p;
            
            final_bitmap = malloc(size_y * size_x / 4);
            if (final_bitmap == NULL) {
                fprintf(stderr, "Error: Unable to allocate memory to generate sprites\n");
                exit(1);
            }
            p = final_bitmap;
            for (c = 0; c < size_y; c += 16) {
                for (d = 0; d < size_x; d += 8) {
                    memcpy(p, bitmap + c / 8 * size_x * 2 + d * 2, 16);
                    p += 16;
                    memcpy(p, bitmap + (c / 8 + 1) * size_x * 2 + d * 2, 16);
                    p += 16;
                }
            }
            fprintf(output, "\t'\n");
            fprintf(output, "%s:\n", label);
            for (c = 0; c < p - final_bitmap; c += 32) {
                generate_data(output, &final_bitmap[c], 8, 32, 0, 0);
                fprintf(output, "\n");
            }
            free(final_bitmap);
        } else if (sms_mode) { /* Sega Master System */
            unsigned char *final_bitmap;
            unsigned char *p;
            
            final_bitmap = malloc(size_y * size_x / 2);
            if (final_bitmap == NULL) {
                fprintf(stderr, "Error: Unable to allocate memory to generate sprites\n");
                exit(1);
            }
            p = final_bitmap;
            for (c = 0; c < size_y; c += 16) {
                for (d = 0; d < size_x; d += 8) {
                    memcpy(p, bitmap + c / 8 * size_x * 4 + d * 4, 32);
                    p += 32;
                    memcpy(p, bitmap + (c / 8 + 1) * size_x * 4 + d * 4, 32);
                    p += 32;
                }
            }
            if (cvbasic) {
                if (sprite_mode == 2) {
                    if (pletter) {
                        fprintf(stderr, "Warning: Sprite mode with BITMAP statements doesn't compress with Pletter\n");
                    }
                    
                    fprintf(output, "\t'\n");
                    fprintf(output, "\t' Sample code:\n");
                    fprintf(output, "\t' DEFINE SPRITE %d,%d,%s\n", 0, size_y / 16 * (size_x / 8), label);
                    fprintf(output, "\t'\n");
                    fprintf(output, "%s:\n", label);
                    for (c = 0; c < p - final_bitmap; c += 64) {
                        generate_data(output, &final_bitmap[c], 8, 64, 0, 0);
                        fprintf(output, "\n");
                    }
                } else {
                    fprintf(output, "\t'\n");
                    fprintf(output, "\t' Sample code:\n");
                    fprintf(output, "\t' DEFINE SPRITE %s%d,%d,%s\n", pletter ? "PLETTER " : "", 0, size_y / 16 * (size_x / 8), label);
                    fprintf(output, "\t'\n");
                    fprintf(output, "%s:\n", label);
                    generate_data(output, final_bitmap, 8, p - final_bitmap, pletter, 1);
                }
            } else {
                fprintf(output, "\t; Total sprites: %d\n", size_y / 16 * (size_x / 8));
                fprintf(output, "%s:\n", label);
                generate_db(output, final_bitmap, 8, p - final_bitmap, pletter);
            }
            free(final_bitmap);
        } else {    /* TMS9918 */
            unsigned char *final_bitmap;
            unsigned char *p;
            
            final_bitmap = malloc((size_y / 16) * (size_x / 8) * 16);
            if (final_bitmap == NULL) {
                fprintf(stderr, "Error: Unable to allocate memory to generate sprites\n");
                exit(1);
            }
            p = final_bitmap;
            for (c = 0; c < size_y; c += 16) {
                for (d = 0; d < size_x; d += 8) {
                    memcpy(p, bitmap + c / 8 * size_x + d, 8);
                    p += 8;
                    memcpy(p, bitmap + (c / 8 + 1) * size_x + d, 8);
                    p += 8;
                }
            }
            if (cvbasic) {
                if (sprite_mode == 2) {
                    if (pletter) {
                        fprintf(stderr, "Warning: Sprite mode with BITMAP statements doesn't compress with Pletter\n");
                    }
                    
                    fprintf(output, "\t'\n");
                    fprintf(output, "\t' Sample code:\n");
                    fprintf(output, "\t' DEFINE SPRITE %d,%d,%s\n", 0, size_y / 16 * (size_x / 16), label);
                    fprintf(output, "\t'\n");
                    fprintf(output, "%s:\n", label);
                    for (c = 0; c < p - final_bitmap; c += 32) {
                        for (d = 0; d < 16; d++) {
                            fprintf(output, "\tBITMAP \"%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c\"\n",
                                    (final_bitmap[c + d] & 0x80) ? 'X' : '.',
                                    (final_bitmap[c + d] & 0x40) ? 'X' : '.',
                                    (final_bitmap[c + d] & 0x20) ? 'X' : '.',
                                    (final_bitmap[c + d] & 0x10) ? 'X' : '.',
                                    (final_bitmap[c + d] & 0x08) ? 'X' : '.',
                                    (final_bitmap[c + d] & 0x04) ? 'X' : '.',
                                    (final_bitmap[c + d] & 0x02) ? 'X' : '.',
                                    (final_bitmap[c + d] & 0x01) ? 'X' : '.',
                                    (final_bitmap[c + d + 16] & 0x80) ? 'X' : '.',
                                    (final_bitmap[c + d + 16] & 0x40) ? 'X' : '.',
                                    (final_bitmap[c + d + 16] & 0x20) ? 'X' : '.',
                                    (final_bitmap[c + d + 16] & 0x10) ? 'X' : '.',
                                    (final_bitmap[c + d + 16] & 0x08) ? 'X' : '.',
                                    (final_bitmap[c + d + 16] & 0x04) ? 'X' : '.',
                                    (final_bitmap[c + d + 16] & 0x02) ? 'X' : '.',
                                    (final_bitmap[c + d + 16] & 0x01) ? 'X' : '.');
                        }
                        fprintf(output, "\n");
                    }
                } else {
                    fprintf(output, "\t'\n");
                    fprintf(output, "\t' Sample code:\n");
                    fprintf(output, "\t' DEFINE SPRITE %s%d,%d,%s\n", pletter ? "PLETTER " : "", 0, size_y / 16 * (size_x / 16), label);
                    fprintf(output, "\t'\n");
                    fprintf(output, "%s:\n", label);
                    generate_data(output, final_bitmap, 8, p - final_bitmap, pletter, 1);
                }
            } else {
                fprintf(output, "\t; Total sprites: %d\n", size_y / 16 * (size_x / 16));
                fprintf(output, "%s:\n", label);
                generate_db(output, final_bitmap, 8, p - final_bitmap, pletter);
            }
            free(final_bitmap);
        }
    } else {    /* Full graphic mode */
        if (cvbasic) {
            if (nes_mode) {
                fprintf(output, "%s_bitmap:\n", label);
                generate_data(output, bitmap, 8, size_x * size_y / 4, pletter, 0);
            } else if (sms_mode) {
                fprintf(output, "%s_bitmap:\n", label);
                generate_data(output, bitmap, 8, size_x * size_y / 2, pletter, 0);
            } else {
                if (size_x * size_y / 8 == 0x1800) {
                    fprintf(output, "\t%sMODE 1\n", remove_stub ? "' " : "");
                    fprintf(output, "\t%sSCREEN DISABLE\n", remove_stub ? "' " : "");
                    fprintf(output, "\t%sDEFINE VRAM %s$0000,$1800,%s_bitmap\n", remove_stub ? "' " : "", pletter ? "PLETTER " : "", label);
                    fprintf(output, "\t%sDEFINE VRAM %s$2000,$1800,%s_color\n", remove_stub ? "' " : "", pletter ? "PLETTER " : "", label);
                    if (magic_sprites)
                        fprintf(output, "\t%sGOSUB %s_show\n", remove_stub ? "' " : "", label);
                    fprintf(output, "\t%sSCREEN ENABLE\n", remove_stub ? "' " : "");
                    if (remove_stub)
                        fprintf(output, "\n");
                    else
                        fprintf(output, "\tWHILE 1: WEND\n\n");
                }
                fprintf(output, "%s_bitmap:\n", label);
                generate_data(output, bitmap, 8, size_x * size_y / 8, pletter, 0);
                fprintf(output, "\n");
                fprintf(output, "%s_color:\n", label);
                generate_data(output, color, 8, size_x * size_y / 8, pletter, 1);
            }
        } else {
            if (nes_mode) {
                fprintf(output, "%s_bitmap:\n", label);
                generate_db(output, bitmap, 8, size_x * size_y / 4, pletter);
            } else if (sms_mode) {
                fprintf(output, "%s_bitmap:\n", label);
                generate_db(output, bitmap, 8, size_x * size_y / 2, pletter);
            } else {
                fprintf(output, "%s_bitmap:\n", label);
                generate_db(output, bitmap, 8, size_x * size_y / 8, pletter);
                fprintf(output, "\n");
                fprintf(output, "%s_color:\n", label);
                generate_db(output, color, 8, size_x * size_y / 8, pletter);
            }
        }
    }
    /* Only works for TMS9918, never executed for Sega Master System */
    if (magic_sprites) {
        if (cvbasic) {
            fprintf(output, "%s_sprites:\n", label);
            generate_data(output, sprites, 8, sig_sprite * 32, pletter, 1);
            fprintf(output, "\n");
            fprintf(output, "%s_show:\tPROCEDURE\n", label);
            if (sig_sprite) {
                fprintf(output, "\tDEFINE SPRITE %s0,%d,%s_sprites\n", pletter ? "PLETTER " : "", sig_sprite, label);
                for (c = 0; c < 128; c += 4) {
                    if (attr[c] != 0xd1)
                        fprintf(output, "\tSPRITE %d,%d,%d,%d,%d\n", c / 4, attr[c], attr[c + 1], attr[c + 2], attr[c + 3]);
                }
            }
            fprintf(output, "\tEND\n");
        } else {
            fprintf(output, "%s_sprites:\n", label);
            generate_db(output, sprites, 8, sizeof(sprites), pletter);
            fprintf(output, "\n");
            fprintf(output, "%s_sprites_attr:\n", label);
            generate_db(output, attr, 4, sizeof(attr), pletter);
        }
    }
    fclose(output);
    free(image);
    return 0;
}
