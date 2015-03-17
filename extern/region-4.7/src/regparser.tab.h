
/* A Bison parser, made by GNU Bison 2.4.1.  */

/* Skeleton interface for Bison's Yacc-like parsers in C
   
      Copyright (C) 1984, 1989, 1990, 2000, 2001, 2002, 2003, 2004, 2005, 2006
   Free Software Foundation, Inc.
   
   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.
   
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
   
   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.  */

/* As a special exception, you may create a larger work that contains
   part or all of the Bison parser skeleton and distribute that work
   under terms of your choice, so long as that work isn't itself a
   parser generator using the skeleton or a modified version thereof
   as a parser skeleton.  Alternatively, if you modify or redistribute
   the parser skeleton itself, you may (at your option) remove this
   special exception, which will cause the skeleton and the resulting
   Bison output files to be licensed under the GNU General Public
   License without this special exception.
   
   This special exception was added by the Free Software Foundation in
   version 2.2 of Bison.  */


/* Tokens.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
   /* Put the tokens into the symbol table, so that GDB and other debuggers
      know about them.  */
   enum yytokentype {
     REG_CIRCLE = 258,
     REG_BOX = 259,
     REG_POLY = 260,
     REG_PIE = 261,
     REG_RECT = 262,
     REG_ELL = 263,
     REG_ROTBOX = 264,
     REG_ANNULUS = 265,
     REG_SECTOR = 266,
     REG_POINT = 267,
     REG_ROTRECT = 268,
     REG_NEG = 269,
     REG_AND = 270,
     REG_OR = 271,
     REG_MINUS = 272,
     REG_FIELD = 273,
     REG_TEXT = 274,
     REG_STRING = 275,
     REG_LINE = 276,
     REG_NUMBER = 277,
     colName = 278
   };
#endif
/* Tokens.  */
#define REG_CIRCLE 258
#define REG_BOX 259
#define REG_POLY 260
#define REG_PIE 261
#define REG_RECT 262
#define REG_ELL 263
#define REG_ROTBOX 264
#define REG_ANNULUS 265
#define REG_SECTOR 266
#define REG_POINT 267
#define REG_ROTRECT 268
#define REG_NEG 269
#define REG_AND 270
#define REG_OR 271
#define REG_MINUS 272
#define REG_FIELD 273
#define REG_TEXT 274
#define REG_STRING 275
#define REG_LINE 276
#define REG_NUMBER 277
#define colName 278




#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
typedef union YYSTYPE
{

/* Line 1676 of yacc.c  */
#line 23 "regparser.y"

  double dval;
  char str[1024];
  regRegion* my_region;
  regShape* my_shape;
  struct polyside
  { double *polyX; double *polyY; long polyS; } PolySide;



/* Line 1676 of yacc.c  */
#line 109 "regparser.tab.h"
} YYSTYPE;
# define YYSTYPE_IS_TRIVIAL 1
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
#endif

extern YYSTYPE yylval;


