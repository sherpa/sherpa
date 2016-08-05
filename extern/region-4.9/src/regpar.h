/*                                                                
**  Copyright (C) 2007  Smithsonian Astrophysical Observatory 
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


                                                
#ifndef reg_PARSER_H
#define reg_PARSER_H
                                             
#ifdef __cplusplus
extern "C" {
#endif


/* The following macros help ensure that the reg parser coexists   */
/* peacefully with other parsers that might be present in the app */
#define yylval    regYYlval


extern int      regYYparse( void );
extern int      regYYlex( void );
void     regYYerror(char* message);
void     regLEXerror(const char* message);
extern int               regYYdebug;  /* Declaration to keep gcc happy */

#define YYDEBUG 1

#ifdef DEBUG
#define YYDEBUG 1
#endif

/*-------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif

#endif










