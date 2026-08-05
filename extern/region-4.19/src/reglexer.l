%{
/*                                                                
**  Copyright (C) 2007,2022  Smithsonian Astrophysical Observatory 
*/                                                                

/*                                                                          */
/*  This program is free software; you can redistribute it and/or modify    */
/*  it under the terms of the GNU General Public License as published by    */
/*  the Free Software Foundation; either version 3 of the License, or       */
/*  (at your option) any later version.                                     */
/*                                                                          */
/*  This program is distributed in the hope that it will be useful,         */
/*  but WITHOUT ANY WARRANTY; without even the implied warranty of          */
/*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           */
/*  GNU General Public License for more details.                            */
/*                                                                          */
/*  You should have received a copy of the GNU General Public License along */
/*  with this program; if not, write to the Free Software Foundation, Inc., */
/*  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.             */
/*                                                                          */

  
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#include "region_priv.h"
/* #include "regionio.h" */
#include "regpar.h"
#include "regparser.h"
char     *regParseStr;
char     *regParseStrEnd;
regRegion *my_Gregion;
void     regYYrestart(FILE*);

  #ifdef FLEX_SCANNER
  #define YY_INPUT(b,r,ms)      (r = reg_flex_input(b,ms))
  #define YY_NO_UNPUT
  int reg_flex_input(char *buf, int max_size);
  #else
  #undef input
  #define input()               (*regParseStr++)
  #define yyless(n)             (regParseStr -= (yyleng - (n)))
  #endif

  /* general lex routines */                 

  #undef  unput                          
  #define unput(c)              (*--regParseStr = c)
  void regLEXerror(const char* message);
%}

%option noyywrap
%option prefix="regYY"

%%

["][0-9A-Za-z~_./#]*["] { strcpy( regYYlval.str, yytext ); return( REG_STRING ); }

[cC][iI][rR]([cC]([Ll]([Ee])?)?)? { return( REG_CIRCLE ); }

[bB][Oo][Xx] { return( REG_BOX ); }
[tT][eE][Xx][Tt] { return( REG_TEXT ); }

[aA][nN][nN][uU][lL][uU][sS] { return ( REG_ANNULUS ); }

[pP][Oo][Ll][Yy]([Gg]([Oo]([Nn])?)?)? { return( REG_POLY); }

[pP][Ii][Ee] { return( REG_PIE ); }

[Ss][Ee][Cc][Tt]([Oo]([Rr])?)? { return( REG_SECTOR ); }

[Rr][Ee][Cc][Tt]([Aa]([Nn]([Gg]([Ll]([Ee])?)?)?)?)? { return( REG_RECT );}

[Ee][Ll][Ll]([Ii]([Pp]([Ss]([Ee])?)?)?)? { return( REG_ELL ); }

[Ll][Ii][Nn][Ee] { return (REG_LINE); }

[Rr][Oo][Tt][Bb][Oo][Xx]  { return( REG_ROTBOX ); }

[Pp][Oo][Ii]([Nn]([Tt])?)? { return( REG_POINT); }

[Rr][Oo][Tt][Rr][Ee][Cc][Tt]([Aa]([Nn]([Gg]([Ll]([Ee])?)?)?)?)? { return( REG_ROTRECT ); }

[Ff][Ii][Ee][Ll][Dd]  { return( REG_FIELD ); }

[Ee][Ll][Ll][Ii][Pp][Tt][Aa][Nn][Nn]([Uu]([Ll]([Uu]([Ss])?)?)?)? { return ( REG_ELLIPTANNULUS ); }

[Dd][Ii][Aa]([Mm]([Oo]([Nn]([Dd])?)?)?)? { return REG_DIAMOND; }

[Rr][Hh][Oo]([Mm]([Bb]([Uu]([Ss])?)?)?)? { return REG_DIAMOND; }

[Rr][Oo][Tt][Dd][Ii][Aa]([Mm]([Oo]([Nn]([Dd])?)?)?)? { return REG_ROTDIAMOND; }

[Rr][Oo][Tt][Rr][Hh][Oo]([Mm]([Bb]([Uu]([Ss])?)?)?)? { return REG_ROTDIAMOND; }

[ ]*[+]?[0-9]+  { regYYlval.dval=atof(yytext); return ( REG_NUMBER ); }
[ ]*[+]?[0-9]+\.?[0-9]*([EeDd][+-]?[0-9]+)?  { regYYlval.dval=atof(yytext); return ( REG_NUMBER ); }
[ ]*[+]?\.[0-9]*([EeDd][+-]?[0-9]+)?  { regYYlval.dval=atof(yytext); return ( REG_NUMBER ); }


[-] { return( REG_MINUS ); }
[!] { return( REG_NEG ); }
[+|] { return( REG_OR ); }
[*&] { return( REG_AND); }

[ \t] ;

\n { return( 0 ); }
.   { return( yytext[0]);} ;

%%


#ifdef FLEX_SCANNER
int
reg_flex_input(char *buf, int max_size)
{
        int n = regParseStrEnd - regParseStr;
        if ( n > max_size ) n = max_size;
        if (n > 0)
        {
                memcpy(buf,regParseStr,n);
                regParseStr += n;
        }
        return n;
}
#endif


