
/* A Bison parser, made by GNU Bison 2.4.1.  */

/* Skeleton implementation for Bison's Yacc-like parsers in C
   
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

/* C LALR(1) parser skeleton written by Richard Stallman, by
   simplifying the original so-called "semantic" parser.  */

/* All symbols defined below should begin with yy or YY, to avoid
   infringing on user name space.  This should be done even for local
   variables, as they might otherwise be expanded by user macros.
   There are some unavoidable exceptions within include files to
   define necessary library symbols; they are noted "INFRINGES ON
   USER NAME SPACE" below.  */

/* Identify Bison output.  */
#define YYBISON 1

/* Bison version.  */
#define YYBISON_VERSION "2.4.1"

/* Skeleton name.  */
#define YYSKELETON_NAME "yacc.c"

/* Pure parsers.  */
#define YYPURE 0

/* Push parsers.  */
#define YYPUSH 0

/* Pull parsers.  */
#define YYPULL 1

/* Using locations.  */
#define YYLSP_NEEDED 0

/* Substitute the variable and function names.  */
#define yyparse         regYYparse
#define yylex           regYYlex
#define yyerror         regYYerror
#define yylval          regYYlval
#define yychar          regYYchar
#define yydebug         regYYdebug
#define yynerrs         regYYnerrs


/* Copy the first part of user declarations.  */

/* Line 189 of yacc.c  */
#line 2 "regparser.y"

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

#include "region_priv.h"
#include "regpar.h"
#include <ctype.h>
#include <float.h>

extern char     *regParseStr;
extern char     *regParseStrEnd;
extern regRegion *my_Gregion;
static int world_coord;
static int world_size;
int test_link_jcm( void );


/* Line 189 of yacc.c  */
#line 116 "regparser.c"

/* Enabling traces.  */
#ifndef YYDEBUG
# define YYDEBUG 0
#endif

/* Enabling verbose error messages.  */
#ifdef YYERROR_VERBOSE
# undef YYERROR_VERBOSE
# define YYERROR_VERBOSE 1
#else
# define YYERROR_VERBOSE 0
#endif

/* Enabling the token table.  */
#ifndef YYTOKEN_TABLE
# define YYTOKEN_TABLE 0
#endif


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
     colName = 278,
     REG_ELLIPTANNULUS = 279,
     REG_DIAMOND = 280,
     REG_ROTDIAMOND = 281
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
#define REG_ELLIPTANNULUS 279
#define REG_DIAMOND 280
#define REG_ROTDIAMOND 281




#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
typedef union YYSTYPE
{

/* Line 214 of yacc.c  */
#line 45 "regparser.y"

  double dval;
  char str[1024];
  regRegion* my_region;
  regShape* my_shape;
  struct polyside
  { double *polyX; double *polyY; long polyS; } PolySide;



/* Line 214 of yacc.c  */
#line 215 "regparser.c"
} YYSTYPE;
# define YYSTYPE_IS_TRIVIAL 1
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
#endif


/* Copy the second part of user declarations.  */


/* Line 264 of yacc.c  */
#line 227 "regparser.c"

#ifdef short
# undef short
#endif

#ifdef YYTYPE_UINT8
typedef YYTYPE_UINT8 yytype_uint8;
#else
typedef unsigned char yytype_uint8;
#endif

#ifdef YYTYPE_INT8
typedef YYTYPE_INT8 yytype_int8;
#elif (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
typedef signed char yytype_int8;
#else
typedef short int yytype_int8;
#endif

#ifdef YYTYPE_UINT16
typedef YYTYPE_UINT16 yytype_uint16;
#else
typedef unsigned short int yytype_uint16;
#endif

#ifdef YYTYPE_INT16
typedef YYTYPE_INT16 yytype_int16;
#else
typedef short int yytype_int16;
#endif

#ifndef YYSIZE_T
# ifdef __SIZE_TYPE__
#  define YYSIZE_T __SIZE_TYPE__
# elif defined size_t
#  define YYSIZE_T size_t
# elif ! defined YYSIZE_T && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
#  include <stddef.h> /* INFRINGES ON USER NAME SPACE */
#  define YYSIZE_T size_t
# else
#  define YYSIZE_T unsigned int
# endif
#endif

#define YYSIZE_MAXIMUM ((YYSIZE_T) -1)

#ifndef YY_
# if YYENABLE_NLS
#  if ENABLE_NLS
#   include <libintl.h> /* INFRINGES ON USER NAME SPACE */
#   define YY_(msgid) dgettext ("bison-runtime", msgid)
#  endif
# endif
# ifndef YY_
#  define YY_(msgid) msgid
# endif
#endif

/* Suppress unused-variable warnings by "using" E.  */
#if ! defined lint || defined __GNUC__
# define YYUSE(e) ((void) (e))
#else
# define YYUSE(e) /* empty */
#endif

/* Identity function, used to suppress warnings about constant conditions.  */
#ifndef lint
# define YYID(n) (n)
#else
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static int
YYID (int yyi)
#else
static int
YYID (yyi)
    int yyi;
#endif
{
  return yyi;
}
#endif

#if ! defined yyoverflow || YYERROR_VERBOSE

/* The parser invokes alloca or malloc; define the necessary symbols.  */

# ifdef YYSTACK_USE_ALLOCA
#  if YYSTACK_USE_ALLOCA
#   ifdef __GNUC__
#    define YYSTACK_ALLOC __builtin_alloca
#   elif defined __BUILTIN_VA_ARG_INCR
#    include <alloca.h> /* INFRINGES ON USER NAME SPACE */
#   elif defined _AIX
#    define YYSTACK_ALLOC __alloca
#   elif defined _MSC_VER
#    include <malloc.h> /* INFRINGES ON USER NAME SPACE */
#    define alloca _alloca
#   else
#    define YYSTACK_ALLOC alloca
#    if ! defined _ALLOCA_H && ! defined _STDLIB_H && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
#     include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#     ifndef _STDLIB_H
#      define _STDLIB_H 1
#     endif
#    endif
#   endif
#  endif
# endif

# ifdef YYSTACK_ALLOC
   /* Pacify GCC's `empty if-body' warning.  */
#  define YYSTACK_FREE(Ptr) do { /* empty */; } while (YYID (0))
#  ifndef YYSTACK_ALLOC_MAXIMUM
    /* The OS might guarantee only one guard page at the bottom of the stack,
       and a page size can be as small as 4096 bytes.  So we cannot safely
       invoke alloca (N) if N exceeds 4096.  Use a slightly smaller number
       to allow for a few compiler-allocated temporary stack slots.  */
#   define YYSTACK_ALLOC_MAXIMUM 4032 /* reasonable circa 2006 */
#  endif
# else
#  define YYSTACK_ALLOC YYMALLOC
#  define YYSTACK_FREE YYFREE
#  ifndef YYSTACK_ALLOC_MAXIMUM
#   define YYSTACK_ALLOC_MAXIMUM YYSIZE_MAXIMUM
#  endif
#  if (defined __cplusplus && ! defined _STDLIB_H \
       && ! ((defined YYMALLOC || defined malloc) \
	     && (defined YYFREE || defined free)))
#   include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#   ifndef _STDLIB_H
#    define _STDLIB_H 1
#   endif
#  endif
#  ifndef YYMALLOC
#   define YYMALLOC malloc
#   if ! defined malloc && ! defined _STDLIB_H && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
void *malloc (YYSIZE_T); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
#  ifndef YYFREE
#   define YYFREE free
#   if ! defined free && ! defined _STDLIB_H && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
void free (void *); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
# endif
#endif /* ! defined yyoverflow || YYERROR_VERBOSE */


#if (! defined yyoverflow \
     && (! defined __cplusplus \
	 || (defined YYSTYPE_IS_TRIVIAL && YYSTYPE_IS_TRIVIAL)))

/* A type that is properly aligned for any stack member.  */
union yyalloc
{
  yytype_int16 yyss_alloc;
  YYSTYPE yyvs_alloc;
};

/* The size of the maximum gap between one aligned stack and the next.  */
# define YYSTACK_GAP_MAXIMUM (sizeof (union yyalloc) - 1)

/* The size of an array large to enough to hold all stacks, each with
   N elements.  */
# define YYSTACK_BYTES(N) \
     ((N) * (sizeof (yytype_int16) + sizeof (YYSTYPE)) \
      + YYSTACK_GAP_MAXIMUM)

/* Copy COUNT objects from FROM to TO.  The source and destination do
   not overlap.  */
# ifndef YYCOPY
#  if defined __GNUC__ && 1 < __GNUC__
#   define YYCOPY(To, From, Count) \
      __builtin_memcpy (To, From, (Count) * sizeof (*(From)))
#  else
#   define YYCOPY(To, From, Count)		\
      do					\
	{					\
	  YYSIZE_T yyi;				\
	  for (yyi = 0; yyi < (Count); yyi++)	\
	    (To)[yyi] = (From)[yyi];		\
	}					\
      while (YYID (0))
#  endif
# endif

/* Relocate STACK from its old location to the new one.  The
   local variables YYSIZE and YYSTACKSIZE give the old and new number of
   elements in the stack, and YYPTR gives the new location of the
   stack.  Advance YYPTR to a properly aligned location for the next
   stack.  */
# define YYSTACK_RELOCATE(Stack_alloc, Stack)				\
    do									\
      {									\
	YYSIZE_T yynewbytes;						\
	YYCOPY (&yyptr->Stack_alloc, Stack, yysize);			\
	Stack = &yyptr->Stack_alloc;					\
	yynewbytes = yystacksize * sizeof (*Stack) + YYSTACK_GAP_MAXIMUM; \
	yyptr += yynewbytes / sizeof (*yyptr);				\
      }									\
    while (YYID (0))

#endif

/* YYFINAL -- State number of the termination state.  */
#define YYFINAL  54
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   588

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  37
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  9
/* YYNRULES -- Number of rules.  */
#define YYNRULES  65
/* YYNRULES -- Number of states.  */
#define YYNSTATES  358

/* YYTRANSLATE(YYLEX) -- Bison symbol number corresponding to YYLEX.  */
#define YYUNDEFTOK  2
#define YYMAXUTOK   282

#define YYTRANSLATE(YYX)						\
  ((unsigned int) (YYX) <= YYMAXUTOK ? yytranslate[YYX] : YYUNDEFTOK)

/* YYTRANSLATE[YYLEX] -- Bison symbol number corresponding to YYLEX.  */
static const yytype_uint8 yytranslate[] =
{
       0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,    31,     2,     2,     2,     2,    30,
      34,    35,     2,     2,    36,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,    29,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
      28,     2,     2,     2,     2,    33,     2,     2,     2,     2,
       2,     2,    32,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     1,     2,     3,     4,
       5,     6,     7,     8,     9,    10,    11,    12,    13,    14,
      15,    16,    17,    18,    19,    20,    21,    22,    23,    24,
      25,    26,    27
};

#if YYDEBUG
/* YYPRHS[YYN] -- Index of the first RHS symbol of rule number YYN in
   YYRHS.  */
static const yytype_uint16 yyprhs[] =
{
       0,     0,     3,     7,    10,    14,    18,    21,    23,    26,
      28,    31,    37,    39,    42,    48,    55,    57,    60,    62,
      65,    68,    71,    73,    76,    79,    82,    85,    89,    92,
      96,   101,   110,   119,   129,   140,   152,   163,   175,   188,
     202,   215,   229,   244,   260,   271,   283,   290,   298,   309,
     321,   332,   344,   357,   371,   376,   382,   397,   414,   430,
     448,   459,   471,   484,   498,   504
};

/* YYRHS -- A `-1'-separated list of the rules' RHS.  */
static const yytype_int8 yyrhs[] =
{
      38,     0,    -1,    38,    16,    44,    -1,    38,    44,    -1,
      38,    15,    44,    -1,    38,    17,    44,    -1,    16,    44,
      -1,    44,    -1,     1,    27,    -1,    41,    -1,    41,    28,
      -1,    22,    29,    22,    29,    22,    -1,    41,    -1,    41,
      28,    -1,    22,    29,    22,    29,    22,    -1,    17,    22,
      29,    22,    29,    22,    -1,    22,    -1,    17,    22,    -1,
      41,    -1,    41,    30,    -1,    41,    31,    -1,    41,    28,
      -1,    41,    -1,    41,    32,    -1,    41,    33,    -1,    41,
      28,    -1,    41,    30,    -1,    41,    30,    30,    -1,    41,
      31,    -1,    18,    34,    35,    -1,    14,    18,    34,    35,
      -1,    19,    34,    39,    36,    40,    36,    20,    35,    -1,
       3,    34,    39,    36,    40,    36,    43,    35,    -1,    14,
       3,    34,    39,    36,    40,    36,    43,    35,    -1,    10,
      34,    39,    36,    40,    36,    43,    36,    43,    35,    -1,
      14,    10,    34,    39,    36,    40,    36,    43,    36,    43,
      35,    -1,     4,    34,    39,    36,    40,    36,    43,    36,
      43,    35,    -1,    14,     4,    34,    39,    36,    40,    36,
      43,    36,    43,    35,    -1,     4,    34,    39,    36,    40,
      36,    43,    36,    43,    36,    42,    35,    -1,    14,     4,
      34,    39,    36,    40,    36,    43,    36,    43,    36,    42,
      35,    -1,     9,    34,    39,    36,    40,    36,    43,    36,
      43,    36,    42,    35,    -1,    14,     9,    34,    39,    36,
      40,    36,    43,    36,    43,    36,    42,    35,    -1,     6,
      34,    39,    36,    40,    36,    43,    36,    43,    36,    42,
      36,    42,    35,    -1,    14,     6,    34,    39,    36,    40,
      36,    43,    36,    43,    36,    42,    36,    42,    35,    -1,
      11,    34,    39,    36,    40,    36,    42,    36,    42,    35,
      -1,    14,    11,    34,    39,    36,    40,    36,    42,    36,
      42,    35,    -1,    12,    34,    39,    36,    40,    35,    -1,
      14,    12,    34,    39,    36,    40,    35,    -1,     7,    34,
      39,    36,    40,    36,    39,    36,    40,    35,    -1,    14,
       7,    34,    39,    36,    40,    36,    39,    36,    40,    35,
      -1,    21,    34,    39,    36,    40,    36,    39,    36,    40,
      35,    -1,    14,    21,    34,    39,    36,    40,    36,    39,
      36,    40,    35,    -1,     8,    34,    39,    36,    40,    36,
      43,    36,    43,    36,    42,    35,    -1,    14,     8,    34,
      39,    36,    40,    36,    43,    36,    43,    36,    42,    35,
      -1,     5,    34,    45,    35,    -1,    14,     5,    34,    45,
      35,    -1,    24,    34,    39,    36,    40,    36,    43,    36,
      43,    36,    43,    36,    43,    35,    -1,    24,    34,    39,
      36,    40,    36,    43,    36,    43,    36,    43,    36,    43,
      36,    42,    35,    -1,    14,    24,    34,    39,    36,    40,
      36,    43,    36,    43,    36,    43,    36,    43,    35,    -1,
      14,    24,    34,    39,    36,    40,    36,    43,    36,    43,
      36,    43,    36,    43,    36,    42,    35,    -1,    25,    34,
      39,    36,    40,    36,    43,    36,    43,    35,    -1,    14,
      25,    34,    39,    36,    40,    36,    43,    36,    43,    35,
      -1,    26,    34,    39,    36,    40,    36,    43,    36,    43,
      36,    42,    35,    -1,    14,    25,    34,    39,    36,    40,
      36,    43,    36,    43,    36,    42,    35,    -1,    45,    36,
      39,    36,    40,    -1,    39,    36,    40,    -1
};

/* YYRLINE[YYN] -- source line where rule number YYN was defined.  */
static const yytype_uint16 yyrline[] =
{
       0,    66,    66,    67,    68,    69,    70,    75,    80,    86,
      88,    90,   101,   103,   105,   113,   124,   126,   131,   133,
     135,   137,   142,   144,   146,   148,   150,   152,   154,   161,
     169,   177,   188,   204,   219,   240,   261,   277,   292,   307,
     322,   337,   352,   371,   390,   400,   410,   420,   430,   440,
     450,   461,   471,   486,   502,   514,   528,   537,   546,   555,
     567,   576,   585,   594,   606,   615
};
#endif

#if YYDEBUG || YYERROR_VERBOSE || YYTOKEN_TABLE
/* YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals.  */
static const char *const yytname[] =
{
  "$end", "error", "$undefined", "REG_CIRCLE", "REG_BOX", "REG_POLY",
  "REG_PIE", "REG_RECT", "REG_ELL", "REG_ROTBOX", "REG_ANNULUS",
  "REG_SECTOR", "REG_POINT", "REG_ROTRECT", "REG_NEG", "REG_AND", "REG_OR",
  "REG_MINUS", "REG_FIELD", "REG_TEXT", "REG_STRING", "REG_LINE",
  "REG_NUMBER", "colName", "REG_ELLIPTANNULUS", "REG_DIAMOND",
  "REG_ROTDIAMOND", "\"\"", "'d'", "':'", "'\\''", "'\"'", "'p'", "'i'",
  "'('", "')'", "','", "$accept", "reg_comp", "reg_xcoord", "reg_ycoord",
  "reg_value", "reg_angle", "reg_size", "reg_shape", "reg_ordered_pair", 0
};
#endif

# ifdef YYPRINT
/* YYTOKNUM[YYLEX-NUM] -- Internal token number corresponding to
   token YYLEX-NUM.  */
static const yytype_uint16 yytoknum[] =
{
       0,   256,   257,   258,   259,   260,   261,   262,   263,   264,
     265,   266,   267,   268,   269,   270,   271,   272,   273,   274,
     275,   276,   277,   278,   279,   280,   281,   282,   100,    58,
      39,    34,   112,   105,    40,    41,    44
};
# endif

/* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const yytype_uint8 yyr1[] =
{
       0,    37,    38,    38,    38,    38,    38,    38,    38,    39,
      39,    39,    40,    40,    40,    40,    41,    41,    42,    42,
      42,    42,    43,    43,    43,    43,    43,    43,    43,    44,
      44,    44,    44,    44,    44,    44,    44,    44,    44,    44,
      44,    44,    44,    44,    44,    44,    44,    44,    44,    44,
      44,    44,    44,    44,    44,    44,    44,    44,    44,    44,
      44,    44,    44,    44,    45,    45
};

/* YYR2[YYN] -- Number of symbols composing right hand side of rule YYN.  */
static const yytype_uint8 yyr2[] =
{
       0,     2,     3,     2,     3,     3,     2,     1,     2,     1,
       2,     5,     1,     2,     5,     6,     1,     2,     1,     2,
       2,     2,     1,     2,     2,     2,     2,     3,     2,     3,
       4,     8,     8,     9,    10,    11,    10,    11,    12,    13,
      12,    13,    14,    15,    10,    11,     6,     7,    10,    11,
      10,    11,    12,    13,     4,     5,    14,    16,    15,    17,
      10,    11,    12,    13,     5,     3
};

/* YYDEFACT[STATE-NAME] -- Default rule to reduce with in state
   STATE-NUM when YYTABLE doesn't specify something else to do.  Zero
   means the default is an error.  */
static const yytype_uint8 yydefact[] =
{
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     7,     8,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     6,     0,     0,
       0,     0,     0,     0,     1,     0,     0,     0,     3,     0,
      16,     0,     9,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,    29,     0,     0,
       0,     0,     0,     4,     2,     5,    17,     0,     0,    10,
       0,     0,    54,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,    30,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,    16,     0,    12,     0,    65,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,    55,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,    17,     0,     0,    13,     0,     0,
       0,     0,     0,     0,     0,     0,    46,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,    11,     0,     0,    16,    22,     0,
       0,    64,     0,     0,     0,     0,     0,    18,     0,     0,
       0,     0,     0,     0,     0,     0,     0,    47,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,    25,    26,
      28,    23,    24,    32,     0,     0,     0,     0,     0,     0,
      21,    19,    20,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,    31,     0,     0,     0,     0,
       0,    14,    27,     0,     0,     0,     0,     0,     0,     0,
      33,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,    15,    36,     0,     0,    48,
       0,     0,    34,    44,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,    50,     0,    60,     0,     0,     0,
       0,     0,    37,     0,     0,    49,     0,     0,    35,    45,
      51,     0,    61,     0,     0,     0,    38,     0,    52,    40,
       0,     0,     0,     0,     0,     0,     0,    62,     0,    39,
       0,    53,    41,     0,    63,     0,    42,     0,     0,    56,
       0,    43,    58,     0,     0,     0,    57,    59
};

/* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int16 yydefgoto[] =
{
      -1,    20,    64,   133,   198,   208,   199,    21,    65
};

/* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
   STATE-NUM.  */
#define YYPACT_NINF -32
static const yytype_int16 yypact[] =
{
     538,    -2,   -13,    -1,    15,    26,    35,    36,    37,    38,
      39,    40,    89,   562,    42,    45,    54,    55,    56,    57,
     505,   -32,   -32,    -6,    -6,    -6,    -6,    -6,    -6,    -6,
      -6,    -6,    -6,    74,    75,    77,    78,    81,    82,    83,
      84,    85,    86,    87,    90,   101,   102,   -32,   103,    -6,
      -6,    -6,    -6,    -6,   -32,   562,   562,   562,   -32,    14,
      10,   104,    18,   105,   106,   -16,   107,   108,   109,   111,
     113,   114,   115,    -6,    -6,    -6,    -6,    -6,    -6,    -6,
      -6,    -6,    -6,   118,    -6,    -6,    -6,   -32,   119,   120,
     121,   122,   123,   -32,   -32,   -32,   -32,   117,     0,   -32,
       0,     0,   -32,    -6,     0,     0,     0,     0,     0,     0,
       0,   124,   125,     2,   126,   127,   128,   129,   130,   132,
     133,   -32,   134,   135,   136,     0,     0,     0,     0,     0,
     146,   161,   155,   151,   158,   152,   -32,   154,   156,   160,
     162,   168,   169,   170,   159,     0,     0,   -32,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   171,   172,
     173,   174,   175,   192,   186,   194,     1,   -32,     1,     0,
       1,    -6,     1,     1,     1,     1,   -32,   181,   182,   183,
     185,   189,   201,   202,   203,   187,   204,   205,   206,   177,
      -6,     1,     1,     1,   -32,   221,   215,   -32,   -18,   210,
     211,   -32,   212,   213,   216,   222,   225,     4,   226,     1,
       1,     1,    -6,     1,     1,     1,     1,   -32,    -6,     1,
       1,   231,   227,   233,   234,   235,   243,   229,   -32,   244,
     -32,   -32,   -32,   -32,     1,     1,     0,     1,     1,     1,
     -32,   -32,   -32,     1,   238,   239,   241,   242,   245,   251,
     252,   253,   259,   261,   262,   -32,     0,     1,     1,     1,
     257,   -32,   -32,     5,   263,   266,   267,   269,   271,   272,
     -32,     1,     1,     0,     1,     1,     1,     1,     0,     1,
       1,   273,   274,   276,   277,   -32,   -32,     1,     1,   -32,
       1,     1,   -32,   -32,     7,   278,   280,   282,   284,   281,
     288,   289,   290,     9,   -32,     1,   -32,     1,   296,   303,
     306,   307,   -32,     1,     1,   -32,     1,     1,   -32,   -32,
     -32,     1,   -32,     1,   308,   310,   -32,     1,   -32,   -32,
     311,   312,   314,   315,   316,   323,     1,   -32,   324,   -32,
       1,   -32,   -32,     1,   -32,    12,   -32,   325,    29,   -32,
       1,   -32,   -32,     1,   326,   328,   -32,   -32
};

/* YYPGOTO[NTERM-NUM].  */
static const yytype_int16 yypgoto[] =
{
     -32,   -32,   150,   228,   -23,   -31,   197,    11,    62
};

/* YYTABLE[YYPACT[STATE-NUM]].  What to do in state STATE-NUM.  If
   positive, shift that token.  If negative, reduce the rule which
   number is the opposite.  If zero, do what YYDEFACT says.
   If YYTABLE_NINF, syntax error.  */
#define YYTABLE_NINF -1
static const yytype_uint16 yytable[] =
{
      62,    62,    62,    62,    62,    62,    62,    62,    62,    62,
     228,    59,   229,   230,   231,   232,    60,   131,    59,   102,
     103,    23,   132,   197,    47,    22,    62,    62,    62,    62,
      62,    58,   240,    24,   241,   242,    96,   147,   103,    97,
     286,   287,   312,   313,   322,   323,    99,   349,   350,    25,
      62,    62,    62,    62,    62,    62,    62,    62,    62,    62,
      26,    62,    62,    62,   352,   353,    93,    94,    95,    27,
      28,    29,    30,    31,    32,   134,    48,   134,   134,    49,
      62,   134,   134,   134,   134,   134,   134,   134,    50,    51,
      52,    53,    33,    34,    35,    36,    37,    38,    39,    40,
      41,    42,   134,   134,   134,   134,   134,    43,    73,    74,
      44,    75,    76,    45,    46,    77,    78,    79,    80,    81,
      82,    83,   134,   134,    84,   134,   134,   134,   134,   134,
     134,   134,   134,   134,   134,    85,    86,   113,    87,   130,
      98,   100,   101,   104,   105,   106,   134,   107,    62,   108,
     109,   110,   207,   121,     0,   125,   126,   127,   128,   129,
     145,   146,   148,   149,   150,   151,   152,    62,   153,   154,
     155,   156,   157,    61,    63,   163,    66,    67,    68,    69,
      70,    71,    72,   164,   165,   251,   167,   166,   168,    62,
     169,     0,   170,   207,   176,    62,   171,   221,   172,    88,
      89,    90,    91,    92,   173,   174,   175,   189,   190,   191,
     192,   193,   269,   134,   194,   195,   196,   209,   210,   211,
     207,   212,   217,   111,   112,   213,   114,   115,   116,   117,
     118,   119,   120,   134,   122,   123,   124,   214,   215,   216,
     218,   219,   220,   226,   227,   233,   300,   234,   235,   236,
     134,   261,   237,   137,   207,   134,   308,   309,   238,   310,
     311,   239,   243,   256,   207,   207,   255,   207,   207,   257,
     258,   259,   260,   270,   262,   271,   325,   272,   273,   285,
       0,   274,   330,   331,   207,   332,   333,   275,   276,   277,
     207,   207,   335,   207,   207,   278,   338,   279,   280,   288,
     207,   289,     0,   290,   207,   291,   292,   293,   304,   347,
     305,   306,     0,   307,   314,   315,   318,   207,   316,   354,
     317,   203,   355,   319,   320,     0,   321,   207,   135,   136,
     207,   326,   138,   139,   140,   141,   142,   143,   144,   327,
     222,   328,   329,     0,   336,   337,   339,     0,   340,   341,
     342,     0,   343,   158,   159,   160,   161,   162,   344,   346,
     351,   356,   247,   357,     0,   200,     0,   202,   252,   204,
     205,   206,     0,   177,   178,     0,   179,   180,   181,   182,
     183,   184,   185,   186,   187,   188,     0,     0,   223,   224,
     225,     0,     0,     0,     0,     0,     0,   201,     0,     0,
       0,     0,     0,     0,     0,     0,   244,   245,   246,     0,
     248,   249,   250,     0,     0,     0,   253,   254,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,   263,   264,     0,   266,   267,   268,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   282,   283,   284,     0,     0,     0,
       0,     0,     0,     0,   265,     0,     0,     0,   294,   295,
       0,   297,   298,   299,     0,     0,   302,   303,     0,     0,
       0,     0,     0,     0,   281,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,   296,   324,     0,     0,    54,   301,     0,     2,     3,
       4,     5,     6,     7,     8,     9,    10,    11,   334,    12,
      55,    56,    57,    14,    15,     0,    16,     0,     0,    17,
      18,    19,     0,   345,     0,     0,     0,     0,     0,     1,
     348,     2,     3,     4,     5,     6,     7,     8,     9,    10,
      11,     0,    12,     0,    13,     0,    14,    15,     0,    16,
       0,     0,    17,    18,    19,     2,     3,     4,     5,     6,
       7,     8,     9,    10,    11,     0,    12,     0,     0,     0,
      14,    15,     0,    16,     0,     0,    17,    18,    19
};

static const yytype_int16 yycheck[] =
{
      23,    24,    25,    26,    27,    28,    29,    30,    31,    32,
      28,    17,    30,    31,    32,    33,    22,    17,    17,    35,
      36,    34,    22,    22,    13,    27,    49,    50,    51,    52,
      53,    20,    28,    34,    30,    31,    22,    35,    36,    29,
      35,    36,    35,    36,    35,    36,    28,    35,    36,    34,
      73,    74,    75,    76,    77,    78,    79,    80,    81,    82,
      34,    84,    85,    86,    35,    36,    55,    56,    57,    34,
      34,    34,    34,    34,    34,    98,    34,   100,   101,    34,
     103,   104,   105,   106,   107,   108,   109,   110,    34,    34,
      34,    34,     3,     4,     5,     6,     7,     8,     9,    10,
      11,    12,   125,   126,   127,   128,   129,    18,    34,    34,
      21,    34,    34,    24,    25,    34,    34,    34,    34,    34,
      34,    34,   145,   146,    34,   148,   149,   150,   151,   152,
     153,   154,   155,   156,   157,    34,    34,    75,    35,    22,
      36,    36,    36,    36,    36,    36,   169,    36,   171,    36,
      36,    36,   175,    35,    -1,    36,    36,    36,    36,    36,
      36,    36,    36,    36,    36,    36,    36,   190,    36,    36,
      36,    36,    36,    23,    24,    29,    26,    27,    28,    29,
      30,    31,    32,    22,    29,   216,    28,    36,    36,   212,
      36,    -1,    36,   216,    35,   218,    36,    20,    36,    49,
      50,    51,    52,    53,    36,    36,    36,    36,    36,    36,
      36,    36,   243,   236,    22,    29,    22,    36,    36,    36,
     243,    36,    35,    73,    74,    36,    76,    77,    78,    79,
      80,    81,    82,   256,    84,    85,    86,    36,    36,    36,
      36,    36,    36,    22,    29,    35,   277,    36,    36,    36,
     273,    22,    36,   103,   277,   278,   287,   288,    36,   290,
     291,    36,    36,    36,   287,   288,    35,   290,   291,    36,
      36,    36,    29,    35,    30,    36,   307,    36,    36,    22,
      -1,    36,   313,   314,   307,   316,   317,    36,    36,    36,
     313,   314,   323,   316,   317,    36,   327,    36,    36,    36,
     323,    35,    -1,    36,   327,    36,    35,    35,    35,   340,
      36,    35,    -1,    36,    36,    35,    35,   340,    36,   350,
      36,   171,   353,    35,    35,    -1,    36,   350,   100,   101,
     353,    35,   104,   105,   106,   107,   108,   109,   110,    36,
     190,    35,    35,    -1,    36,    35,    35,    -1,    36,    35,
      35,    -1,    36,   125,   126,   127,   128,   129,    35,    35,
      35,    35,   212,    35,    -1,   168,    -1,   170,   218,   172,
     173,   174,    -1,   145,   146,    -1,   148,   149,   150,   151,
     152,   153,   154,   155,   156,   157,    -1,    -1,   191,   192,
     193,    -1,    -1,    -1,    -1,    -1,    -1,   169,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,   209,   210,   211,    -1,
     213,   214,   215,    -1,    -1,    -1,   219,   220,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,   234,   235,    -1,   237,   238,   239,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,   257,   258,   259,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,   236,    -1,    -1,    -1,   271,   272,
      -1,   274,   275,   276,    -1,    -1,   279,   280,    -1,    -1,
      -1,    -1,    -1,    -1,   256,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,   273,   305,    -1,    -1,     0,   278,    -1,     3,     4,
       5,     6,     7,     8,     9,    10,    11,    12,   321,    14,
      15,    16,    17,    18,    19,    -1,    21,    -1,    -1,    24,
      25,    26,    -1,   336,    -1,    -1,    -1,    -1,    -1,     1,
     343,     3,     4,     5,     6,     7,     8,     9,    10,    11,
      12,    -1,    14,    -1,    16,    -1,    18,    19,    -1,    21,
      -1,    -1,    24,    25,    26,     3,     4,     5,     6,     7,
       8,     9,    10,    11,    12,    -1,    14,    -1,    -1,    -1,
      18,    19,    -1,    21,    -1,    -1,    24,    25,    26
};

/* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
   symbol of state STATE-NUM.  */
static const yytype_uint8 yystos[] =
{
       0,     1,     3,     4,     5,     6,     7,     8,     9,    10,
      11,    12,    14,    16,    18,    19,    21,    24,    25,    26,
      38,    44,    27,    34,    34,    34,    34,    34,    34,    34,
      34,    34,    34,     3,     4,     5,     6,     7,     8,     9,
      10,    11,    12,    18,    21,    24,    25,    44,    34,    34,
      34,    34,    34,    34,     0,    15,    16,    17,    44,    17,
      22,    39,    41,    39,    39,    45,    39,    39,    39,    39,
      39,    39,    39,    34,    34,    34,    34,    34,    34,    34,
      34,    34,    34,    34,    34,    34,    34,    35,    39,    39,
      39,    39,    39,    44,    44,    44,    22,    29,    36,    28,
      36,    36,    35,    36,    36,    36,    36,    36,    36,    36,
      36,    39,    39,    45,    39,    39,    39,    39,    39,    39,
      39,    35,    39,    39,    39,    36,    36,    36,    36,    36,
      22,    17,    22,    40,    41,    40,    40,    39,    40,    40,
      40,    40,    40,    40,    40,    36,    36,    35,    36,    36,
      36,    36,    36,    36,    36,    36,    36,    36,    40,    40,
      40,    40,    40,    29,    22,    29,    36,    28,    36,    36,
      36,    36,    36,    36,    36,    36,    35,    40,    40,    40,
      40,    40,    40,    40,    40,    40,    40,    40,    40,    36,
      36,    36,    36,    36,    22,    29,    22,    22,    41,    43,
      43,    40,    43,    39,    43,    43,    43,    41,    42,    36,
      36,    36,    36,    36,    36,    36,    36,    35,    36,    36,
      36,    20,    39,    43,    43,    43,    22,    29,    28,    30,
      31,    32,    33,    35,    36,    36,    36,    36,    36,    36,
      28,    30,    31,    36,    43,    43,    43,    39,    43,    43,
      43,    42,    39,    43,    43,    35,    36,    36,    36,    36,
      29,    22,    30,    43,    43,    40,    43,    43,    43,    42,
      35,    36,    36,    36,    36,    36,    36,    36,    36,    36,
      36,    40,    43,    43,    43,    22,    35,    36,    36,    35,
      36,    36,    35,    35,    43,    43,    40,    43,    43,    43,
      42,    40,    43,    43,    35,    36,    35,    36,    42,    42,
      42,    42,    35,    36,    36,    35,    36,    36,    35,    35,
      35,    36,    35,    36,    43,    42,    35,    36,    35,    35,
      42,    42,    42,    42,    43,    42,    36,    35,    42,    35,
      36,    35,    35,    36,    35,    43,    35,    42,    43,    35,
      36,    35,    35,    36,    42,    42,    35,    35
};

#define yyerrok		(yyerrstatus = 0)
#define yyclearin	(yychar = YYEMPTY)
#define YYEMPTY		(-2)
#define YYEOF		0

#define YYACCEPT	goto yyacceptlab
#define YYABORT		goto yyabortlab
#define YYERROR		goto yyerrorlab


/* Like YYERROR except do call yyerror.  This remains here temporarily
   to ease the transition to the new meaning of YYERROR, for GCC.
   Once GCC version 2 has supplanted version 1, this can go.  */

#define YYFAIL		goto yyerrlab

#define YYRECOVERING()  (!!yyerrstatus)

#define YYBACKUP(Token, Value)					\
do								\
  if (yychar == YYEMPTY && yylen == 1)				\
    {								\
      yychar = (Token);						\
      yylval = (Value);						\
      yytoken = YYTRANSLATE (yychar);				\
      YYPOPSTACK (1);						\
      goto yybackup;						\
    }								\
  else								\
    {								\
      yyerror (YY_("syntax error: cannot back up")); \
      YYERROR;							\
    }								\
while (YYID (0))


#define YYTERROR	1
#define YYERRCODE	256


/* YYLLOC_DEFAULT -- Set CURRENT to span from RHS[1] to RHS[N].
   If N is 0, then set CURRENT to the empty location which ends
   the previous symbol: RHS[0] (always defined).  */

#define YYRHSLOC(Rhs, K) ((Rhs)[K])
#ifndef YYLLOC_DEFAULT
# define YYLLOC_DEFAULT(Current, Rhs, N)				\
    do									\
      if (YYID (N))                                                    \
	{								\
	  (Current).first_line   = YYRHSLOC (Rhs, 1).first_line;	\
	  (Current).first_column = YYRHSLOC (Rhs, 1).first_column;	\
	  (Current).last_line    = YYRHSLOC (Rhs, N).last_line;		\
	  (Current).last_column  = YYRHSLOC (Rhs, N).last_column;	\
	}								\
      else								\
	{								\
	  (Current).first_line   = (Current).last_line   =		\
	    YYRHSLOC (Rhs, 0).last_line;				\
	  (Current).first_column = (Current).last_column =		\
	    YYRHSLOC (Rhs, 0).last_column;				\
	}								\
    while (YYID (0))
#endif


/* YY_LOCATION_PRINT -- Print the location on the stream.
   This macro was not mandated originally: define only if we know
   we won't break user code: when these are the locations we know.  */

#ifndef YY_LOCATION_PRINT
# if YYLTYPE_IS_TRIVIAL
#  define YY_LOCATION_PRINT(File, Loc)			\
     fprintf (File, "%d.%d-%d.%d",			\
	      (Loc).first_line, (Loc).first_column,	\
	      (Loc).last_line,  (Loc).last_column)
# else
#  define YY_LOCATION_PRINT(File, Loc) ((void) 0)
# endif
#endif


/* YYLEX -- calling `yylex' with the right arguments.  */

#ifdef YYLEX_PARAM
# define YYLEX yylex (YYLEX_PARAM)
#else
# define YYLEX yylex ()
#endif

/* Enable debugging if requested.  */
#if YYDEBUG

# ifndef YYFPRINTF
#  include <stdio.h> /* INFRINGES ON USER NAME SPACE */
#  define YYFPRINTF fprintf
# endif

# define YYDPRINTF(Args)			\
do {						\
  if (yydebug)					\
    YYFPRINTF Args;				\
} while (YYID (0))

# define YY_SYMBOL_PRINT(Title, Type, Value, Location)			  \
do {									  \
  if (yydebug)								  \
    {									  \
      YYFPRINTF (stderr, "%s ", Title);					  \
      yy_symbol_print (stderr,						  \
		  Type, Value); \
      YYFPRINTF (stderr, "\n");						  \
    }									  \
} while (YYID (0))


/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
`--------------------------------*/

/*ARGSUSED*/
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_symbol_value_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep)
#else
static void
yy_symbol_value_print (yyoutput, yytype, yyvaluep)
    FILE *yyoutput;
    int yytype;
    YYSTYPE const * const yyvaluep;
#endif
{
  if (!yyvaluep)
    return;
# ifdef YYPRINT
  if (yytype < YYNTOKENS)
    YYPRINT (yyoutput, yytoknum[yytype], *yyvaluep);
# else
  YYUSE (yyoutput);
# endif
  switch (yytype)
    {
      default:
	break;
    }
}


/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
`--------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_symbol_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep)
#else
static void
yy_symbol_print (yyoutput, yytype, yyvaluep)
    FILE *yyoutput;
    int yytype;
    YYSTYPE const * const yyvaluep;
#endif
{
  if (yytype < YYNTOKENS)
    YYFPRINTF (yyoutput, "token %s (", yytname[yytype]);
  else
    YYFPRINTF (yyoutput, "nterm %s (", yytname[yytype]);

  yy_symbol_value_print (yyoutput, yytype, yyvaluep);
  YYFPRINTF (yyoutput, ")");
}

/*------------------------------------------------------------------.
| yy_stack_print -- Print the state stack from its BOTTOM up to its |
| TOP (included).                                                   |
`------------------------------------------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_stack_print (yytype_int16 *yybottom, yytype_int16 *yytop)
#else
static void
yy_stack_print (yybottom, yytop)
    yytype_int16 *yybottom;
    yytype_int16 *yytop;
#endif
{
  YYFPRINTF (stderr, "Stack now");
  for (; yybottom <= yytop; yybottom++)
    {
      int yybot = *yybottom;
      YYFPRINTF (stderr, " %d", yybot);
    }
  YYFPRINTF (stderr, "\n");
}

# define YY_STACK_PRINT(Bottom, Top)				\
do {								\
  if (yydebug)							\
    yy_stack_print ((Bottom), (Top));				\
} while (YYID (0))


/*------------------------------------------------.
| Report that the YYRULE is going to be reduced.  |
`------------------------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_reduce_print (YYSTYPE *yyvsp, int yyrule)
#else
static void
yy_reduce_print (yyvsp, yyrule)
    YYSTYPE *yyvsp;
    int yyrule;
#endif
{
  int yynrhs = yyr2[yyrule];
  int yyi;
  unsigned long int yylno = yyrline[yyrule];
  YYFPRINTF (stderr, "Reducing stack by rule %d (line %lu):\n",
	     yyrule - 1, yylno);
  /* The symbols being reduced.  */
  for (yyi = 0; yyi < yynrhs; yyi++)
    {
      YYFPRINTF (stderr, "   $%d = ", yyi + 1);
      yy_symbol_print (stderr, yyrhs[yyprhs[yyrule] + yyi],
		       &(yyvsp[(yyi + 1) - (yynrhs)])
		       		       );
      YYFPRINTF (stderr, "\n");
    }
}

# define YY_REDUCE_PRINT(Rule)		\
do {					\
  if (yydebug)				\
    yy_reduce_print (yyvsp, Rule); \
} while (YYID (0))

/* Nonzero means print parse trace.  It is left uninitialized so that
   multiple parsers can coexist.  */
int yydebug;
#else /* !YYDEBUG */
# define YYDPRINTF(Args)
# define YY_SYMBOL_PRINT(Title, Type, Value, Location)
# define YY_STACK_PRINT(Bottom, Top)
# define YY_REDUCE_PRINT(Rule)
#endif /* !YYDEBUG */


/* YYINITDEPTH -- initial size of the parser's stacks.  */
#ifndef	YYINITDEPTH
# define YYINITDEPTH 200
#endif

/* YYMAXDEPTH -- maximum size the stacks can grow to (effective only
   if the built-in stack extension method is used).

   Do not make this value too large; the results are undefined if
   YYSTACK_ALLOC_MAXIMUM < YYSTACK_BYTES (YYMAXDEPTH)
   evaluated with infinite-precision integer arithmetic.  */

#ifndef YYMAXDEPTH
# define YYMAXDEPTH 10000
#endif



#if YYERROR_VERBOSE

# ifndef yystrlen
#  if defined __GLIBC__ && defined _STRING_H
#   define yystrlen strlen
#  else
/* Return the length of YYSTR.  */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static YYSIZE_T
yystrlen (const char *yystr)
#else
static YYSIZE_T
yystrlen (yystr)
    const char *yystr;
#endif
{
  YYSIZE_T yylen;
  for (yylen = 0; yystr[yylen]; yylen++)
    continue;
  return yylen;
}
#  endif
# endif

# ifndef yystpcpy
#  if defined __GLIBC__ && defined _STRING_H && defined _GNU_SOURCE
#   define yystpcpy stpcpy
#  else
/* Copy YYSRC to YYDEST, returning the address of the terminating '\0' in
   YYDEST.  */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static char *
yystpcpy (char *yydest, const char *yysrc)
#else
static char *
yystpcpy (yydest, yysrc)
    char *yydest;
    const char *yysrc;
#endif
{
  char *yyd = yydest;
  const char *yys = yysrc;

  while ((*yyd++ = *yys++) != '\0')
    continue;

  return yyd - 1;
}
#  endif
# endif

# ifndef yytnamerr
/* Copy to YYRES the contents of YYSTR after stripping away unnecessary
   quotes and backslashes, so that it's suitable for yyerror.  The
   heuristic is that double-quoting is unnecessary unless the string
   contains an apostrophe, a comma, or backslash (other than
   backslash-backslash).  YYSTR is taken from yytname.  If YYRES is
   null, do not copy; instead, return the length of what the result
   would have been.  */
static YYSIZE_T
yytnamerr (char *yyres, const char *yystr)
{
  if (*yystr == '"')
    {
      YYSIZE_T yyn = 0;
      char const *yyp = yystr;

      for (;;)
	switch (*++yyp)
	  {
	  case '\'':
	  case ',':
	    goto do_not_strip_quotes;

	  case '\\':
	    if (*++yyp != '\\')
	      goto do_not_strip_quotes;
	    /* Fall through.  */
	  default:
	    if (yyres)
	      yyres[yyn] = *yyp;
	    yyn++;
	    break;

	  case '"':
	    if (yyres)
	      yyres[yyn] = '\0';
	    return yyn;
	  }
    do_not_strip_quotes: ;
    }

  if (! yyres)
    return yystrlen (yystr);

  return yystpcpy (yyres, yystr) - yyres;
}
# endif

/* Copy into YYRESULT an error message about the unexpected token
   YYCHAR while in state YYSTATE.  Return the number of bytes copied,
   including the terminating null byte.  If YYRESULT is null, do not
   copy anything; just return the number of bytes that would be
   copied.  As a special case, return 0 if an ordinary "syntax error"
   message will do.  Return YYSIZE_MAXIMUM if overflow occurs during
   size calculation.  */
static YYSIZE_T
yysyntax_error (char *yyresult, int yystate, int yychar)
{
  int yyn = yypact[yystate];

  if (! (YYPACT_NINF < yyn && yyn <= YYLAST))
    return 0;
  else
    {
      int yytype = YYTRANSLATE (yychar);
      YYSIZE_T yysize0 = yytnamerr (0, yytname[yytype]);
      YYSIZE_T yysize = yysize0;
      YYSIZE_T yysize1;
      int yysize_overflow = 0;
      enum { YYERROR_VERBOSE_ARGS_MAXIMUM = 5 };
      char const *yyarg[YYERROR_VERBOSE_ARGS_MAXIMUM];
      int yyx;

# if 0
      /* This is so xgettext sees the translatable formats that are
	 constructed on the fly.  */
      YY_("syntax error, unexpected %s");
      YY_("syntax error, unexpected %s, expecting %s");
      YY_("syntax error, unexpected %s, expecting %s or %s");
      YY_("syntax error, unexpected %s, expecting %s or %s or %s");
      YY_("syntax error, unexpected %s, expecting %s or %s or %s or %s");
# endif
      char *yyfmt;
      char const *yyf;
      static char const yyunexpected[] = "syntax error, unexpected %s";
      static char const yyexpecting[] = ", expecting %s";
      static char const yyor[] = " or %s";
      char yyformat[sizeof yyunexpected
		    + sizeof yyexpecting - 1
		    + ((YYERROR_VERBOSE_ARGS_MAXIMUM - 2)
		       * (sizeof yyor - 1))];
      char const *yyprefix = yyexpecting;

      /* Start YYX at -YYN if negative to avoid negative indexes in
	 YYCHECK.  */
      int yyxbegin = yyn < 0 ? -yyn : 0;

      /* Stay within bounds of both yycheck and yytname.  */
      int yychecklim = YYLAST - yyn + 1;
      int yyxend = yychecklim < YYNTOKENS ? yychecklim : YYNTOKENS;
      int yycount = 1;

      yyarg[0] = yytname[yytype];
      yyfmt = yystpcpy (yyformat, yyunexpected);

      for (yyx = yyxbegin; yyx < yyxend; ++yyx)
	if (yycheck[yyx + yyn] == yyx && yyx != YYTERROR)
	  {
	    if (yycount == YYERROR_VERBOSE_ARGS_MAXIMUM)
	      {
		yycount = 1;
		yysize = yysize0;
		yyformat[sizeof yyunexpected - 1] = '\0';
		break;
	      }
	    yyarg[yycount++] = yytname[yyx];
	    yysize1 = yysize + yytnamerr (0, yytname[yyx]);
	    yysize_overflow |= (yysize1 < yysize);
	    yysize = yysize1;
	    yyfmt = yystpcpy (yyfmt, yyprefix);
	    yyprefix = yyor;
	  }

      yyf = YY_(yyformat);
      yysize1 = yysize + yystrlen (yyf);
      yysize_overflow |= (yysize1 < yysize);
      yysize = yysize1;

      if (yysize_overflow)
	return YYSIZE_MAXIMUM;

      if (yyresult)
	{
	  /* Avoid sprintf, as that infringes on the user's name space.
	     Don't have undefined behavior even if the translation
	     produced a string with the wrong number of "%s"s.  */
	  char *yyp = yyresult;
	  int yyi = 0;
	  while ((*yyp = *yyf) != '\0')
	    {
	      if (*yyp == '%' && yyf[1] == 's' && yyi < yycount)
		{
		  yyp += yytnamerr (yyp, yyarg[yyi++]);
		  yyf += 2;
		}
	      else
		{
		  yyp++;
		  yyf++;
		}
	    }
	}
      return yysize;
    }
}
#endif /* YYERROR_VERBOSE */


/*-----------------------------------------------.
| Release the memory associated to this symbol.  |
`-----------------------------------------------*/

/*ARGSUSED*/
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yydestruct (const char *yymsg, int yytype, YYSTYPE *yyvaluep)
#else
static void
yydestruct (yymsg, yytype, yyvaluep)
    const char *yymsg;
    int yytype;
    YYSTYPE *yyvaluep;
#endif
{
  YYUSE (yyvaluep);

  if (!yymsg)
    yymsg = "Deleting";
  YY_SYMBOL_PRINT (yymsg, yytype, yyvaluep, yylocationp);

  switch (yytype)
    {

      default:
	break;
    }
}

/* Prevent warnings from -Wmissing-prototypes.  */
#ifdef YYPARSE_PARAM
#if defined __STDC__ || defined __cplusplus
int yyparse (void *YYPARSE_PARAM);
#else
int yyparse ();
#endif
#else /* ! YYPARSE_PARAM */
#if defined __STDC__ || defined __cplusplus
int yyparse (void);
#else
int yyparse ();
#endif
#endif /* ! YYPARSE_PARAM */


/* The lookahead symbol.  */
int yychar;

/* The semantic value of the lookahead symbol.  */
YYSTYPE yylval;

/* Number of syntax errors so far.  */
int yynerrs;



/*-------------------------.
| yyparse or yypush_parse.  |
`-------------------------*/

#ifdef YYPARSE_PARAM
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
int
yyparse (void *YYPARSE_PARAM)
#else
int
yyparse (YYPARSE_PARAM)
    void *YYPARSE_PARAM;
#endif
#else /* ! YYPARSE_PARAM */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
int
yyparse (void)
#else
int
yyparse ()

#endif
#endif
{


    int yystate;
    /* Number of tokens to shift before error messages enabled.  */
    int yyerrstatus;

    /* The stacks and their tools:
       `yyss': related to states.
       `yyvs': related to semantic values.

       Refer to the stacks thru separate pointers, to allow yyoverflow
       to reallocate them elsewhere.  */

    /* The state stack.  */
    yytype_int16 yyssa[YYINITDEPTH];
    yytype_int16 *yyss;
    yytype_int16 *yyssp;

    /* The semantic value stack.  */
    YYSTYPE yyvsa[YYINITDEPTH];
    YYSTYPE *yyvs;
    YYSTYPE *yyvsp;

    YYSIZE_T yystacksize;

  int yyn;
  int yyresult;
  /* Lookahead token as an internal (translated) token number.  */
  int yytoken;
  /* The variables used to return semantic value and location from the
     action routines.  */
  YYSTYPE yyval;

#if YYERROR_VERBOSE
  /* Buffer for error messages, and its allocated size.  */
  char yymsgbuf[128];
  char *yymsg = yymsgbuf;
  YYSIZE_T yymsg_alloc = sizeof yymsgbuf;
#endif

#define YYPOPSTACK(N)   (yyvsp -= (N), yyssp -= (N))

  /* The number of symbols on the RHS of the reduced rule.
     Keep to zero when no symbol should be popped.  */
  int yylen = 0;

  yytoken = 0;
  yyss = yyssa;
  yyvs = yyvsa;
  yystacksize = YYINITDEPTH;

  YYDPRINTF ((stderr, "Starting parse\n"));

  yystate = 0;
  yyerrstatus = 0;
  yynerrs = 0;
  yychar = YYEMPTY; /* Cause a token to be read.  */

  /* Initialize stack pointers.
     Waste one element of value and location stack
     so that they stay on the same level as the state stack.
     The wasted elements are never initialized.  */
  yyssp = yyss;
  yyvsp = yyvs;

  goto yysetstate;

/*------------------------------------------------------------.
| yynewstate -- Push a new state, which is found in yystate.  |
`------------------------------------------------------------*/
 yynewstate:
  /* In all cases, when you get here, the value and location stacks
     have just been pushed.  So pushing a state here evens the stacks.  */
  yyssp++;

 yysetstate:
  *yyssp = yystate;

  if (yyss + yystacksize - 1 <= yyssp)
    {
      /* Get the current used size of the three stacks, in elements.  */
      YYSIZE_T yysize = yyssp - yyss + 1;

#ifdef yyoverflow
      {
	/* Give user a chance to reallocate the stack.  Use copies of
	   these so that the &'s don't force the real ones into
	   memory.  */
	YYSTYPE *yyvs1 = yyvs;
	yytype_int16 *yyss1 = yyss;

	/* Each stack pointer address is followed by the size of the
	   data in use in that stack, in bytes.  This used to be a
	   conditional around just the two extra args, but that might
	   be undefined if yyoverflow is a macro.  */
	yyoverflow (YY_("memory exhausted"),
		    &yyss1, yysize * sizeof (*yyssp),
		    &yyvs1, yysize * sizeof (*yyvsp),
		    &yystacksize);

	yyss = yyss1;
	yyvs = yyvs1;
      }
#else /* no yyoverflow */
# ifndef YYSTACK_RELOCATE
      goto yyexhaustedlab;
# else
      /* Extend the stack our own way.  */
      if (YYMAXDEPTH <= yystacksize)
	goto yyexhaustedlab;
      yystacksize *= 2;
      if (YYMAXDEPTH < yystacksize)
	yystacksize = YYMAXDEPTH;

      {
	yytype_int16 *yyss1 = yyss;
	union yyalloc *yyptr =
	  (union yyalloc *) YYSTACK_ALLOC (YYSTACK_BYTES (yystacksize));
	if (! yyptr)
	  goto yyexhaustedlab;
	YYSTACK_RELOCATE (yyss_alloc, yyss);
	YYSTACK_RELOCATE (yyvs_alloc, yyvs);
#  undef YYSTACK_RELOCATE
	if (yyss1 != yyssa)
	  YYSTACK_FREE (yyss1);
      }
# endif
#endif /* no yyoverflow */

      yyssp = yyss + yysize - 1;
      yyvsp = yyvs + yysize - 1;

      YYDPRINTF ((stderr, "Stack size increased to %lu\n",
		  (unsigned long int) yystacksize));

      if (yyss + yystacksize - 1 <= yyssp)
	YYABORT;
    }

  YYDPRINTF ((stderr, "Entering state %d\n", yystate));

  if (yystate == YYFINAL)
    YYACCEPT;

  goto yybackup;

/*-----------.
| yybackup.  |
`-----------*/
yybackup:

  /* Do appropriate processing given the current state.  Read a
     lookahead token if we need one and don't already have one.  */

  /* First try to decide what to do without reference to lookahead token.  */
  yyn = yypact[yystate];
  if (yyn == YYPACT_NINF)
    goto yydefault;

  /* Not known => get a lookahead token if don't already have one.  */

  /* YYCHAR is either YYEMPTY or YYEOF or a valid lookahead symbol.  */
  if (yychar == YYEMPTY)
    {
      YYDPRINTF ((stderr, "Reading a token: "));
      yychar = YYLEX;
    }

  if (yychar <= YYEOF)
    {
      yychar = yytoken = YYEOF;
      YYDPRINTF ((stderr, "Now at end of input.\n"));
    }
  else
    {
      yytoken = YYTRANSLATE (yychar);
      YY_SYMBOL_PRINT ("Next token is", yytoken, &yylval, &yylloc);
    }

  /* If the proper action on seeing token YYTOKEN is to reduce or to
     detect an error, take that action.  */
  yyn += yytoken;
  if (yyn < 0 || YYLAST < yyn || yycheck[yyn] != yytoken)
    goto yydefault;
  yyn = yytable[yyn];
  if (yyn <= 0)
    {
      if (yyn == 0 || yyn == YYTABLE_NINF)
	goto yyerrlab;
      yyn = -yyn;
      goto yyreduce;
    }

  /* Count tokens shifted since error; after three, turn off error
     status.  */
  if (yyerrstatus)
    yyerrstatus--;

  /* Shift the lookahead token.  */
  YY_SYMBOL_PRINT ("Shifting", yytoken, &yylval, &yylloc);

  /* Discard the shifted token.  */
  yychar = YYEMPTY;

  yystate = yyn;
  *++yyvsp = yylval;

  goto yynewstate;


/*-----------------------------------------------------------.
| yydefault -- do the default action for the current state.  |
`-----------------------------------------------------------*/
yydefault:
  yyn = yydefact[yystate];
  if (yyn == 0)
    goto yyerrlab;
  goto yyreduce;


/*-----------------------------.
| yyreduce -- Do a reduction.  |
`-----------------------------*/
yyreduce:
  /* yyn is the number of a rule to reduce with.  */
  yylen = yyr2[yyn];

  /* If YYLEN is nonzero, implement the default value of the action:
     `$$ = $1'.

     Otherwise, the following line sets YYVAL to garbage.
     This behavior is undocumented and Bison
     users should not rely upon it.  Assigning to YYVAL
     unconditionally makes the parser a bit smaller, and it avoids a
     GCC warning that YYVAL may be used uninitialized.  */
  yyval = yyvsp[1-yylen];


  YY_REDUCE_PRINT (yyn);
  switch (yyn)
    {
        case 2:

/* Line 1455 of yacc.c  */
#line 66 "regparser.y"
    { (yyval.my_region) = (yyvsp[(1) - (3)].my_region); regAddShape( (yyvsp[(1) - (3)].my_region), regOR, (yyvsp[(3) - (3)].my_shape) );  }
    break;

  case 3:

/* Line 1455 of yacc.c  */
#line 67 "regparser.y"
    { (yyval.my_region) = (yyvsp[(1) - (2)].my_region); regAddShape( (yyvsp[(1) - (2)].my_region), regOR, (yyvsp[(2) - (2)].my_shape) );  }
    break;

  case 4:

/* Line 1455 of yacc.c  */
#line 68 "regparser.y"
    { (yyval.my_region) = (yyvsp[(1) - (3)].my_region); regAddShape( (yyvsp[(1) - (3)].my_region), regAND, (yyvsp[(3) - (3)].my_shape) ); }
    break;

  case 5:

/* Line 1455 of yacc.c  */
#line 69 "regparser.y"
    { (yyval.my_region) = (yyvsp[(1) - (3)].my_region); regNegate( (yyvsp[(3) - (3)].my_shape) ); regAddShape( (yyvsp[(1) - (3)].my_region), regAND, (yyvsp[(3) - (3)].my_shape) ); }
    break;

  case 6:

/* Line 1455 of yacc.c  */
#line 70 "regparser.y"
    {
                (yyval.my_region) = regCreateRegion( NULL, NULL );
                regAddShape( (yyval.my_region) , regOR, (yyvsp[(2) - (2)].my_shape) ); 
		my_Gregion = (yyval.my_region);
             }
    break;

  case 7:

/* Line 1455 of yacc.c  */
#line 75 "regparser.y"
    {
                (yyval.my_region) = regCreateRegion(NULL, NULL );
                regAddShape( (yyval.my_region) , regOR, (yyvsp[(1) - (1)].my_shape) ); 
		        my_Gregion = (yyval.my_region);
             }
    break;

  case 8:

/* Line 1455 of yacc.c  */
#line 80 "regparser.y"
    { yyerrok; }
    break;

  case 9:

/* Line 1455 of yacc.c  */
#line 87 "regparser.y"
    {  (yyval.dval) = (yyvsp[(1) - (1)].dval); world_coord = RC_UNK; }
    break;

  case 10:

/* Line 1455 of yacc.c  */
#line 89 "regparser.y"
    {  (yyval.dval) = (yyvsp[(1) - (2)].dval); world_coord = RC_WORLD; }
    break;

  case 11:

/* Line 1455 of yacc.c  */
#line 91 "regparser.y"
    {
           double ndeg, nmin, nsec, nval;
           ndeg = (yyvsp[(1) - (5)].dval); nmin = (yyvsp[(3) - (5)].dval); nsec = (yyvsp[(5) - (5)].dval);
           nval = ( ndeg + nmin / 60.0 + nsec / 3600.0 );
           world_coord = RC_WORLD;
           (yyval.dval) = nval * 15.0;
          }
    break;

  case 12:

/* Line 1455 of yacc.c  */
#line 102 "regparser.y"
    {  (yyval.dval) = (yyvsp[(1) - (1)].dval); world_coord = RC_UNK; }
    break;

  case 13:

/* Line 1455 of yacc.c  */
#line 104 "regparser.y"
    {  (yyval.dval) = (yyvsp[(1) - (2)].dval); world_coord = RC_WORLD; }
    break;

  case 14:

/* Line 1455 of yacc.c  */
#line 106 "regparser.y"
    {
           double ndeg, nmin, nsec, nval;
           ndeg = (yyvsp[(1) - (5)].dval); nmin = (yyvsp[(3) - (5)].dval); nsec = (yyvsp[(5) - (5)].dval);
           nval = ( ndeg + nmin / 60.0 + nsec / 3600.0 );
           world_coord = RC_WORLD;
           (yyval.dval) = nval;
          }
    break;

  case 15:

/* Line 1455 of yacc.c  */
#line 114 "regparser.y"
    {  /* Special care in case of -00:00:01 */
           double ndeg, nmin, nsec, nval;
           ndeg = (yyvsp[(2) - (6)].dval); nmin = (yyvsp[(4) - (6)].dval); nsec = (yyvsp[(6) - (6)].dval);
           nval = ( ndeg + nmin / 60.0 + nsec / 3600.0 );
           world_coord = RC_WORLD;
           (yyval.dval) = -nval;
          }
    break;

  case 16:

/* Line 1455 of yacc.c  */
#line 125 "regparser.y"
    {  (yyval.dval) = (yyvsp[(1) - (1)].dval); }
    break;

  case 17:

/* Line 1455 of yacc.c  */
#line 127 "regparser.y"
    {  (yyval.dval) = -(yyvsp[(2) - (2)].dval);}
    break;

  case 18:

/* Line 1455 of yacc.c  */
#line 132 "regparser.y"
    {  (yyval.dval) = (yyvsp[(1) - (1)].dval); }
    break;

  case 19:

/* Line 1455 of yacc.c  */
#line 134 "regparser.y"
    {  (yyval.dval) = (yyvsp[(1) - (2)].dval) / 60.0; }
    break;

  case 20:

/* Line 1455 of yacc.c  */
#line 136 "regparser.y"
    {  (yyval.dval) = (yyvsp[(1) - (2)].dval) / 3600.0; }
    break;

  case 21:

/* Line 1455 of yacc.c  */
#line 138 "regparser.y"
    {  (yyval.dval) = (yyvsp[(1) - (2)].dval); }
    break;

  case 22:

/* Line 1455 of yacc.c  */
#line 143 "regparser.y"
    {  (yyval.dval) = (yyvsp[(1) - (1)].dval); world_size = RC_UNK; }
    break;

  case 23:

/* Line 1455 of yacc.c  */
#line 145 "regparser.y"
    {  (yyval.dval) = (yyvsp[(1) - (2)].dval); world_size = RC_PHYSICAL; }
    break;

  case 24:

/* Line 1455 of yacc.c  */
#line 147 "regparser.y"
    {  (yyval.dval) = (yyvsp[(1) - (2)].dval); world_size = RC_LOGICAL; }
    break;

  case 25:

/* Line 1455 of yacc.c  */
#line 149 "regparser.y"
    {  (yyval.dval) = (yyvsp[(1) - (2)].dval); world_size = RC_WORLD; }
    break;

  case 26:

/* Line 1455 of yacc.c  */
#line 151 "regparser.y"
    {  (yyval.dval) = (yyvsp[(1) - (2)].dval) / 60.0; world_size = RC_WORLD; }
    break;

  case 27:

/* Line 1455 of yacc.c  */
#line 153 "regparser.y"
    {  (yyval.dval) = (yyvsp[(1) - (3)].dval) / 3600.0; world_size = RC_WORLD; }
    break;

  case 28:

/* Line 1455 of yacc.c  */
#line 155 "regparser.y"
    {  (yyval.dval) = (yyvsp[(1) - (2)].dval) / 3600.0; world_size = RC_WORLD; }
    break;

  case 29:

/* Line 1455 of yacc.c  */
#line 162 "regparser.y"
    { 
       (yyval.my_shape) = regCreateNewWorldShape( regFIELD, regInclude, NULL, NULL, 0, NULL, NULL, world_coord, world_size );
       if ( (yyval.my_shape) == NULL ) {
	 my_Gregion = NULL;
	 YYERROR;
       }
     }
    break;

  case 30:

/* Line 1455 of yacc.c  */
#line 170 "regparser.y"
    { 
       (yyval.my_shape) = regCreateNewWorldShape( regFIELD, regExclude, NULL, NULL, 0, NULL, NULL, world_coord, world_size );
       if ( (yyval.my_shape) == NULL ) {
	 my_Gregion = NULL;
	 YYERROR;
       }
     }
    break;

  case 31:

/* Line 1455 of yacc.c  */
#line 178 "regparser.y"
    { 
       /* We just ignore the text, treat it as a point */
       double x[1]; double y[1];
       x[0]=(yyvsp[(3) - (8)].dval); y[0]=(yyvsp[(5) - (8)].dval);
       (yyval.my_shape) = regCreateNewWorldShape( regPOINT, regInclude, x, y, 1, NULL, NULL, world_coord, 0 );
       if ( (yyval.my_shape) == NULL ) {
	 my_Gregion = NULL;
	 YYERROR;
       }
     }
    break;

  case 32:

/* Line 1455 of yacc.c  */
#line 189 "regparser.y"
    { 
       double x[1]; double y[1]; double r[1];
       x[0]=(yyvsp[(3) - (8)].dval); y[0]=(yyvsp[(5) - (8)].dval); r[0]=(yyvsp[(7) - (8)].dval);

       if ( r[0] < 0 ) {
	 fprintf( stderr, "ERROR: circle radius must be positive\n");
	 my_Gregion = NULL;
	 YYERROR;
       } else
	 (yyval.my_shape) = regCreateNewWorldShape( regCIRCLE, regInclude, x, y, 1, r, NULL, world_coord, world_size );
         if ( (yyval.my_shape) == NULL ) {
	   my_Gregion = NULL;
	   YYERROR;
         }
     }
    break;

  case 33:

/* Line 1455 of yacc.c  */
#line 205 "regparser.y"
    { 
       double x[1]; double y[1]; double r[1];
       x[0]=(yyvsp[(4) - (9)].dval); y[0]=(yyvsp[(6) - (9)].dval); r[0]=(yyvsp[(8) - (9)].dval);
       if ( r[0] < 0 ) {
	 fprintf( stderr, "ERROR: circle radius must be positive\n");
	 my_Gregion = NULL;
	 YYERROR;
       } else
	 (yyval.my_shape) = regCreateNewWorldShape( regCIRCLE, regExclude, x, y, 1, r, NULL, world_coord, world_size );
         if ( (yyval.my_shape) == NULL ) {
  	   my_Gregion = NULL;
	   YYERROR;
         }
     }
    break;

  case 34:

/* Line 1455 of yacc.c  */
#line 220 "regparser.y"
    {
       double x[1]; double y[1]; double r[2];
       x[0]=(yyvsp[(3) - (10)].dval); y[0]=(yyvsp[(5) - (10)].dval); r[0]=(yyvsp[(7) - (10)].dval); r[1]=(yyvsp[(9) - (10)].dval);

       if (( r[0] < 0 ) || (r[1] < 0 )) {
	 fprintf( stderr, "ERROR: annulus radius must be positive\n");
	 my_Gregion = NULL;
	 YYERROR;
       } else if (r[0] > r[1]){
	 fprintf( stderr, "ERROR: outer annuls radius must be >= inner radius\n");
	 my_Gregion = NULL;
	 YYERROR;
       }
       else
	 (yyval.my_shape) = regCreateNewWorldShape( regANNULUS, regInclude, x, y, 1, r, NULL, world_coord, world_size);
         if ( (yyval.my_shape) == NULL ) {
 	   my_Gregion = NULL;
	   YYERROR;
         }
     }
    break;

  case 35:

/* Line 1455 of yacc.c  */
#line 241 "regparser.y"
    {
       double x[1]; double y[1]; double r[2];
       x[0]=(yyvsp[(4) - (11)].dval); y[0]=(yyvsp[(6) - (11)].dval); r[0]=(yyvsp[(8) - (11)].dval); r[1]=(yyvsp[(10) - (11)].dval);

       if (( r[0] < 0 ) || (r[1] < 0 )) {
	 fprintf( stderr, "ERROR: annulus radius must be positive\n");
	 my_Gregion = NULL;
	 YYERROR;
       } else if (r[0] > r[1]){
	 fprintf( stderr, "ERROR: outer annuls radius must be >= inner radius\n");
	 my_Gregion = NULL;
	 YYERROR;
       }
       else
	 (yyval.my_shape) = regCreateNewWorldShape( regANNULUS, regExclude, x, y, 1, r, NULL, world_coord, world_size);
         if ( (yyval.my_shape) == NULL ) {
 	   my_Gregion = NULL;
	   YYERROR;
         }
     }
    break;

  case 36:

/* Line 1455 of yacc.c  */
#line 262 "regparser.y"
    {
       double x[1]; double y[1]; double r[2];
       x[0]=(yyvsp[(3) - (10)].dval); y[0]=(yyvsp[(5) - (10)].dval); r[0]=(yyvsp[(7) - (10)].dval); r[1]=(yyvsp[(9) - (10)].dval);

       if (( r[0] < 0 ) || (r[1] < 0 )) {
	 fprintf( stderr, "ERROR: box lengths must be positive\n");
	 my_Gregion = NULL;
	 YYERROR;
       } else
	 (yyval.my_shape) = regCreateNewWorldShape( regBOX, regInclude, x, y, 1, r, NULL, world_coord, world_size);
         if ( (yyval.my_shape) == NULL ) {
 	   my_Gregion = NULL;
	   YYERROR;
         }
     }
    break;

  case 37:

/* Line 1455 of yacc.c  */
#line 278 "regparser.y"
    {
       double x[1]; double y[1]; double r[2];
       x[0]=(yyvsp[(4) - (11)].dval); y[0]=(yyvsp[(6) - (11)].dval); r[0]=(yyvsp[(8) - (11)].dval); r[1]=(yyvsp[(10) - (11)].dval);
       if (( r[0] < 0 ) || (r[1] < 0 )) {
	 fprintf( stderr, "ERROR: box lengths must be positive\n");
	 my_Gregion = NULL;
	 YYERROR;
       } else
	 (yyval.my_shape) = regCreateNewWorldShape( regBOX, regExclude, x, y, 1, r, NULL, world_coord, world_size);
         if ( (yyval.my_shape) == NULL ) {
 	   my_Gregion = NULL;
	   YYERROR;
         }
     }
    break;

  case 38:

/* Line 1455 of yacc.c  */
#line 293 "regparser.y"
    {
       double x[1]; double y[1]; double r[2]; double a[1];
       x[0]=(yyvsp[(3) - (12)].dval); y[0]=(yyvsp[(5) - (12)].dval); r[0]=(yyvsp[(7) - (12)].dval); r[1]=(yyvsp[(9) - (12)].dval); a[0]=(yyvsp[(11) - (12)].dval);
       if (( r[0] < 0 ) || (r[1] < 0 )) {
	 fprintf( stderr, "ERROR: box lengths must be positive\n");
	 my_Gregion = NULL;
	 YYERROR;
       } else
	 (yyval.my_shape) = regCreateNewWorldShape( regROTBOX, regInclude, x, y, 1, r, a, world_coord, world_size);
         if ( (yyval.my_shape) == NULL ) {
 	   my_Gregion = NULL;
	   YYERROR;
         }
     }
    break;

  case 39:

/* Line 1455 of yacc.c  */
#line 308 "regparser.y"
    {
       double x[1]; double y[1]; double r[2]; double a[1];
       x[0]=(yyvsp[(4) - (13)].dval); y[0]=(yyvsp[(6) - (13)].dval); r[0]=(yyvsp[(8) - (13)].dval); r[1]=(yyvsp[(10) - (13)].dval); a[0]=(yyvsp[(12) - (13)].dval);
       if (( r[0] < 0 ) || (r[1] < 0 )) {
	 fprintf( stderr, "ERROR: box lengths must be positive\n");
	 my_Gregion = NULL;
	 YYERROR;
       } else
	 (yyval.my_shape) = regCreateNewWorldShape( regROTBOX, regExclude, x, y, 1, r, a, world_coord, world_size);
         if ( (yyval.my_shape) == NULL ) {
 	   my_Gregion = NULL;
	   YYERROR;
         }
     }
    break;

  case 40:

/* Line 1455 of yacc.c  */
#line 323 "regparser.y"
    {
       double x[1]; double y[1]; double r[2]; double a[1];
       x[0]=(yyvsp[(3) - (12)].dval); y[0]=(yyvsp[(5) - (12)].dval); r[0]=(yyvsp[(7) - (12)].dval); r[1]=(yyvsp[(9) - (12)].dval); a[0]=(yyvsp[(11) - (12)].dval);
       if (( r[0] < 0 ) || (r[1] < 0 )) {
	 fprintf( stderr, "ERROR: box lengths must be positive\n");
	 my_Gregion = NULL;
	 YYERROR;
       } else
	 (yyval.my_shape) = regCreateNewWorldShape( regROTBOX, regInclude, x, y, 1, r, a, world_coord, world_size);
         if ( (yyval.my_shape) == NULL ) {
 	   my_Gregion = NULL;
	   YYERROR;
         }
     }
    break;

  case 41:

/* Line 1455 of yacc.c  */
#line 338 "regparser.y"
    {
       double x[1]; double y[1]; double r[2]; double a[1];
       x[0]=(yyvsp[(4) - (13)].dval); y[0]=(yyvsp[(6) - (13)].dval); r[0]=(yyvsp[(8) - (13)].dval); r[1]=(yyvsp[(10) - (13)].dval); a[0]=(yyvsp[(12) - (13)].dval);
       if (( r[0] < 0 ) || (r[1] < 0 )) {
	 fprintf( stderr, "ERROR: box lengths must be positive\n");
	 my_Gregion = NULL;
	 YYERROR;
       } else
	 (yyval.my_shape) = regCreateNewWorldShape( regROTBOX, regExclude, x, y, 1, r, a, world_coord, world_size);
         if ( (yyval.my_shape) == NULL ) {
 	   my_Gregion = NULL;
	   YYERROR;
         }
     }
    break;

  case 42:

/* Line 1455 of yacc.c  */
#line 353 "regparser.y"
    {
       double x[1]; double y[1]; double a[2]; double r[2];
       x[0]=(yyvsp[(3) - (14)].dval); y[0]=(yyvsp[(5) - (14)].dval); r[0]=(yyvsp[(7) - (14)].dval); r[1]=(yyvsp[(9) - (14)].dval);  a[0]=(yyvsp[(11) - (14)].dval); a[1]=(yyvsp[(13) - (14)].dval);
       if (( r[0] < 0 ) || (r[1] < 0 )) {
	 fprintf( stderr, "ERROR: pie radii must be positive\n");
	 my_Gregion = NULL;
	 YYERROR;
       } else if ( r[0] > r[1] ) {
	 fprintf( stderr, "ERROR: outer annuls radius must be >= inner radius\n");
	 my_Gregion = NULL;
	 YYERROR;
       } else
	 (yyval.my_shape) = regCreateNewWorldShape( regPIE, regInclude, x, y, 1, r, a, world_coord, world_size );
         if ( (yyval.my_shape) == NULL ) {
 	   my_Gregion = NULL;
	   YYERROR;
         }
     }
    break;

  case 43:

/* Line 1455 of yacc.c  */
#line 372 "regparser.y"
    {
       double x[1]; double y[1]; double a[2]; double r[2];
       x[0]=(yyvsp[(4) - (15)].dval); y[0]=(yyvsp[(6) - (15)].dval); r[0] = (yyvsp[(8) - (15)].dval); r[1] = (yyvsp[(10) - (15)].dval); a[0]=(yyvsp[(12) - (15)].dval); a[1]=(yyvsp[(14) - (15)].dval);
       if (( r[0] < 0 ) || (r[1] < 0 )) {
	 fprintf( stderr, "ERROR: pie radii must be positive\n");
	 my_Gregion = NULL;
	 YYERROR;
       } else if ( r[0] > r[1] ) {
	 fprintf( stderr, "ERROR: outer annuls radius must be >= inner radius\n");
	 my_Gregion = NULL;
	 YYERROR;
       } else
	 (yyval.my_shape) = regCreateNewWorldShape( regPIE, regExclude, x, y, 1, r, a, world_coord, world_size );
         if ( (yyval.my_shape) == NULL ) {
 	   my_Gregion = NULL;
	   YYERROR;
         }
     }
    break;

  case 44:

/* Line 1455 of yacc.c  */
#line 391 "regparser.y"
    {
       double x[1]; double y[1]; double a[2];
       x[0]=(yyvsp[(3) - (10)].dval); y[0]=(yyvsp[(5) - (10)].dval); a[0]=(yyvsp[(7) - (10)].dval); a[1]=(yyvsp[(9) - (10)].dval);
       (yyval.my_shape) = regCreateNewWorldShape( regSECTOR, regInclude, x, y, 1, NULL, a, world_coord, 0 );
       if ( (yyval.my_shape) == NULL ) {
 	 my_Gregion = NULL;
	 YYERROR;
       }
     }
    break;

  case 45:

/* Line 1455 of yacc.c  */
#line 401 "regparser.y"
    {
       double x[1]; double y[1]; double a[2];
       x[0]=(yyvsp[(4) - (11)].dval); y[0]=(yyvsp[(6) - (11)].dval); a[0]=(yyvsp[(8) - (11)].dval); a[1]=(yyvsp[(10) - (11)].dval);
       (yyval.my_shape) = regCreateNewWorldShape( regSECTOR, regExclude, x, y, 1, NULL, a, world_coord, 0 );
       if ( (yyval.my_shape) == NULL ) {
 	 my_Gregion = NULL;
	 YYERROR;
       }
     }
    break;

  case 46:

/* Line 1455 of yacc.c  */
#line 411 "regparser.y"
    {
       double x[1]; double y[1];
       x[0]=(yyvsp[(3) - (6)].dval); y[0]=(yyvsp[(5) - (6)].dval);
       (yyval.my_shape) = regCreateNewWorldShape( regPOINT, regInclude, x, y, 1, NULL, NULL, world_coord, 0 );
       if ( (yyval.my_shape) == NULL ) {
 	 my_Gregion = NULL;
	 YYERROR;
       }
     }
    break;

  case 47:

/* Line 1455 of yacc.c  */
#line 421 "regparser.y"
    {
       double x[1]; double y[1];
       x[0]=(yyvsp[(4) - (7)].dval); y[0]=(yyvsp[(6) - (7)].dval);
       (yyval.my_shape) = regCreateNewWorldShape( regPOINT, regExclude, x, y, 1, NULL, NULL, world_coord, 0 );
       if ( (yyval.my_shape) == NULL ) {
 	 my_Gregion = NULL;
	 YYERROR;
       }
     }
    break;

  case 48:

/* Line 1455 of yacc.c  */
#line 431 "regparser.y"
    {
       double x[2]; double y[2];
       x[0]=(yyvsp[(3) - (10)].dval); x[1]=(yyvsp[(7) - (10)].dval); y[0]=(yyvsp[(5) - (10)].dval);y[1]=(yyvsp[(9) - (10)].dval);
       (yyval.my_shape) = regCreateNewWorldShape( regRECTANGLE, regInclude, x, y, 2, NULL, NULL, world_coord, 0 );
       if ( (yyval.my_shape) == NULL ) {
 	 my_Gregion = NULL;
	 YYERROR;
       }
     }
    break;

  case 49:

/* Line 1455 of yacc.c  */
#line 441 "regparser.y"
    {
       double x[2]; double y[2];
       x[0]=(yyvsp[(4) - (11)].dval); x[1]=(yyvsp[(8) - (11)].dval); y[0]=(yyvsp[(6) - (11)].dval); y[1]=(yyvsp[(10) - (11)].dval);
       (yyval.my_shape) = regCreateNewWorldShape( regRECTANGLE, regExclude, x, y, 2, NULL, NULL, world_coord, 0 );
       if ( (yyval.my_shape) == NULL ) {
 	 my_Gregion = NULL;
	 YYERROR;
       }
     }
    break;

  case 50:

/* Line 1455 of yacc.c  */
#line 451 "regparser.y"
    {
       /* RegLine doesn't work correctly; need to calculate angle */
       double x[2]; double y[2];
       x[0]=(yyvsp[(3) - (10)].dval); x[1]=(yyvsp[(7) - (10)].dval); y[0]=(yyvsp[(5) - (10)].dval);y[1]=(yyvsp[(9) - (10)].dval);
       (yyval.my_shape) = regCreateNewWorldShape( regRECTANGLE, regInclude, x, y, 2, NULL, NULL, world_coord, 0 );
       if ( (yyval.my_shape) == NULL ) {
 	 my_Gregion = NULL;
	 YYERROR;
       }
     }
    break;

  case 51:

/* Line 1455 of yacc.c  */
#line 462 "regparser.y"
    {
       double x[2]; double y[2];
       x[0]=(yyvsp[(4) - (11)].dval); x[1]=(yyvsp[(8) - (11)].dval); y[0]=(yyvsp[(6) - (11)].dval); y[1]=(yyvsp[(10) - (11)].dval);
       (yyval.my_shape) = regCreateNewWorldShape( regRECTANGLE, regExclude, x, y, 2, NULL, NULL, world_coord, 0 );
       if ( (yyval.my_shape) == NULL ) {
 	 my_Gregion = NULL;
	 YYERROR;
       }
     }
    break;

  case 52:

/* Line 1455 of yacc.c  */
#line 472 "regparser.y"
    {
       double x[1]; double y[1]; double r[2]; double a[1];
       x[0]=(yyvsp[(3) - (12)].dval); y[0]=(yyvsp[(5) - (12)].dval); r[0]=(yyvsp[(7) - (12)].dval); r[1]=(yyvsp[(9) - (12)].dval); a[0]=(yyvsp[(11) - (12)].dval);
       if (( r[0] < 0 ) || (r[1] < 0 )) {
	 fprintf( stderr, "ERROR: ellipse radii must be positive\n");
	 my_Gregion = NULL;
	 YYERROR;
       } else
	 (yyval.my_shape) = regCreateNewWorldShape( regELLIPSE, regInclude, x, y, 1, r, a, world_coord, world_size );
         if ( (yyval.my_shape) == NULL ) {
 	   my_Gregion = NULL;
	   YYERROR;
         }
     }
    break;

  case 53:

/* Line 1455 of yacc.c  */
#line 487 "regparser.y"
    {
       double x[1]; double y[1]; double r[2]; double a[1];
       x[0]=(yyvsp[(4) - (13)].dval); y[0]=(yyvsp[(6) - (13)].dval); r[0]=(yyvsp[(8) - (13)].dval); r[1]=(yyvsp[(10) - (13)].dval); a[0]=(yyvsp[(12) - (13)].dval);
       if (( r[0] < 0 ) || (r[1] < 0 )) {
	 fprintf( stderr, "ERROR: ellipse radii must be positive\n");
	 my_Gregion = NULL;
	 YYERROR;
       } else
	 (yyval.my_shape) = regCreateNewWorldShape( regELLIPSE, regExclude, x, y, 1, r, a, world_coord, world_size );
         if ( (yyval.my_shape) == NULL ) {
 	   my_Gregion = NULL;
	   YYERROR;
         }
     }
    break;

  case 54:

/* Line 1455 of yacc.c  */
#line 503 "regparser.y"
    { 
       (yyval.my_shape) = regCreateNewWorldShape( regPOLYGON, regInclude, (yyvsp[(3) - (4)].PolySide).polyX, (yyvsp[(3) - (4)].PolySide).polyY, 
			    (yyvsp[(3) - (4)].PolySide).polyS, NULL, NULL, world_coord, 0 );
       if ( (yyval.my_shape) == NULL ) {
	 my_Gregion = NULL;
	 YYERROR;
       }

       free( (yyvsp[(3) - (4)].PolySide).polyX );
       free( (yyvsp[(3) - (4)].PolySide).polyY );
     }
    break;

  case 55:

/* Line 1455 of yacc.c  */
#line 515 "regparser.y"
    { 
       (yyval.my_shape) = regCreateNewWorldShape( regPOLYGON, regExclude, (yyvsp[(4) - (5)].PolySide).polyX, (yyvsp[(4) - (5)].PolySide).polyY, 
			    (yyvsp[(4) - (5)].PolySide).polyS, NULL, NULL, world_coord, 0 );
       free( (yyvsp[(4) - (5)].PolySide).polyX );
       free( (yyvsp[(4) - (5)].PolySide).polyY );
       if ( (yyval.my_shape) == NULL ) {
	 my_Gregion = NULL;
	 YYERROR;
       }
     }
    break;

  case 56:

/* Line 1455 of yacc.c  */
#line 529 "regparser.y"
    {
        double x[1]; double y[1]; double r1[2]; double r2[2];
        x[0]=(yyvsp[(3) - (14)].dval); y[0]=(yyvsp[(5) - (14)].dval); r1[0]=(yyvsp[(7) - (14)].dval); r1[1]=(yyvsp[(9) - (14)].dval); r2[0] = (yyvsp[(11) - (14)].dval); r2[1] = (yyvsp[(13) - (14)].dval);
 	    
        fprintf(stderr, "ERROR: Elliptannuli are not yet supported.\n");
        my_Gregion = NULL;
	    YYERROR;
     }
    break;

  case 57:

/* Line 1455 of yacc.c  */
#line 538 "regparser.y"
    {
        double x[1]; double y[1]; double r1[2]; double r2[2];  double a[1];
        x[0]=(yyvsp[(3) - (16)].dval); y[0]=(yyvsp[(5) - (16)].dval); r1[0]=(yyvsp[(7) - (16)].dval); r1[1]=(yyvsp[(9) - (16)].dval); r2[0] = (yyvsp[(11) - (16)].dval); r2[1] = (yyvsp[(13) - (16)].dval); a[0]=(yyvsp[(15) - (16)].dval);
 	    
        fprintf(stderr, "ERROR: Elliptannuli are not yet supported.\n");
        my_Gregion = NULL;
	    YYERROR;
     }
    break;

  case 58:

/* Line 1455 of yacc.c  */
#line 547 "regparser.y"
    {
        double x[1]; double y[1]; double r1[2]; double r2[2];
        x[0]=(yyvsp[(4) - (15)].dval); y[0]=(yyvsp[(6) - (15)].dval); r1[0]=(yyvsp[(8) - (15)].dval); r1[1]=(yyvsp[(10) - (15)].dval); r2[0] = (yyvsp[(12) - (15)].dval); r2[1] = (yyvsp[(14) - (15)].dval);
 	    
        fprintf(stderr, "ERROR: Elliptannuli are not yet supported.\n");
        my_Gregion = NULL;
	    YYERROR;
     }
    break;

  case 59:

/* Line 1455 of yacc.c  */
#line 556 "regparser.y"
    {
        double x[1]; double y[1]; double r1[2]; double r2[2];  double a[1];
        x[0]=(yyvsp[(4) - (17)].dval); y[0]=(yyvsp[(6) - (17)].dval); r1[0]=(yyvsp[(8) - (17)].dval); r1[1]=(yyvsp[(10) - (17)].dval); r2[0] = (yyvsp[(12) - (17)].dval); r2[1] = (yyvsp[(14) - (17)].dval); a[0]=(yyvsp[(16) - (17)].dval);
 	    
        fprintf(stderr, "ERROR: Elliptannuli are not yet supported.\n");
        my_Gregion = NULL;
	    YYERROR;
     }
    break;

  case 60:

/* Line 1455 of yacc.c  */
#line 568 "regparser.y"
    {
        double x[1]; double y[1]; double r[2]; 
        x[0]=(yyvsp[(3) - (10)].dval); y[0]=(yyvsp[(5) - (10)].dval); r[0]=(yyvsp[(7) - (10)].dval); r[1]=(yyvsp[(9) - (10)].dval); 
 	    
        fprintf(stderr, "ERROR: Diamonds are not yet supported.\n");
        my_Gregion = NULL;
	    YYERROR;
     }
    break;

  case 61:

/* Line 1455 of yacc.c  */
#line 577 "regparser.y"
    {
        double x[1]; double y[1]; double r[2]; 
        x[0]=(yyvsp[(4) - (11)].dval); y[0]=(yyvsp[(6) - (11)].dval); r[0]=(yyvsp[(8) - (11)].dval); r[1]=(yyvsp[(10) - (11)].dval); 
 	    
        fprintf(stderr, "ERROR: Diamonds are not yet supported.\n");
        my_Gregion = NULL;
	    YYERROR;
     }
    break;

  case 62:

/* Line 1455 of yacc.c  */
#line 586 "regparser.y"
    {
        double x[1]; double y[1]; double r[2]; double a[1];
        x[0]=(yyvsp[(3) - (12)].dval); y[0]=(yyvsp[(5) - (12)].dval); r[0]=(yyvsp[(7) - (12)].dval); r[1]=(yyvsp[(9) - (12)].dval); a[0] = (yyvsp[(11) - (12)].dval);
 	    
        fprintf(stderr, "ERROR: Diamonds are not yet supported.\n");
        my_Gregion = NULL;
	    YYERROR;
     }
    break;

  case 63:

/* Line 1455 of yacc.c  */
#line 595 "regparser.y"
    {
        double x[1]; double y[1]; double r[2]; double a[1];
        x[0]=(yyvsp[(4) - (13)].dval); y[0]=(yyvsp[(6) - (13)].dval); r[0]=(yyvsp[(8) - (13)].dval); r[1]=(yyvsp[(10) - (13)].dval); a[0] = (yyvsp[(12) - (13)].dval);
 	    
        fprintf(stderr, "ERROR: Diamonds are not yet supported.\n");
        my_Gregion = NULL;
	    YYERROR;
     }
    break;

  case 64:

/* Line 1455 of yacc.c  */
#line 607 "regparser.y"
    {
     (yyval.PolySide) = (yyvsp[(1) - (5)].PolySide);
     (yyval.PolySide).polyS += 1;
     (yyval.PolySide).polyX = (double*)realloc( (yyval.PolySide).polyX, (yyval.PolySide).polyS *sizeof(double));
     (yyval.PolySide).polyY = (double*)realloc( (yyval.PolySide).polyY, (yyval.PolySide).polyS *sizeof(double));
     (yyval.PolySide).polyX[(yyval.PolySide).polyS -1] = (yyvsp[(3) - (5)].dval);
     (yyval.PolySide).polyY[(yyval.PolySide).polyS -1] = (yyvsp[(5) - (5)].dval);
     }
    break;

  case 65:

/* Line 1455 of yacc.c  */
#line 615 "regparser.y"
    {
     (yyval.PolySide).polyS = 1;
     (yyval.PolySide).polyX = (double*)calloc(1,sizeof(double));
     (yyval.PolySide).polyY = (double*)calloc(1,sizeof(double));
     (yyval.PolySide).polyX[0] = (yyvsp[(1) - (3)].dval);
     (yyval.PolySide).polyY[0] = (yyvsp[(3) - (3)].dval);
   }
    break;



/* Line 1455 of yacc.c  */
#line 2563 "regparser.c"
      default: break;
    }
  YY_SYMBOL_PRINT ("-> $$ =", yyr1[yyn], &yyval, &yyloc);

  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);

  *++yyvsp = yyval;

  /* Now `shift' the result of the reduction.  Determine what state
     that goes to, based on the state we popped back to and the rule
     number reduced by.  */

  yyn = yyr1[yyn];

  yystate = yypgoto[yyn - YYNTOKENS] + *yyssp;
  if (0 <= yystate && yystate <= YYLAST && yycheck[yystate] == *yyssp)
    yystate = yytable[yystate];
  else
    yystate = yydefgoto[yyn - YYNTOKENS];

  goto yynewstate;


/*------------------------------------.
| yyerrlab -- here on detecting error |
`------------------------------------*/
yyerrlab:
  /* If not already recovering from an error, report this error.  */
  if (!yyerrstatus)
    {
      ++yynerrs;
#if ! YYERROR_VERBOSE
      yyerror (YY_("syntax error"));
#else
      {
	YYSIZE_T yysize = yysyntax_error (0, yystate, yychar);
	if (yymsg_alloc < yysize && yymsg_alloc < YYSTACK_ALLOC_MAXIMUM)
	  {
	    YYSIZE_T yyalloc = 2 * yysize;
	    if (! (yysize <= yyalloc && yyalloc <= YYSTACK_ALLOC_MAXIMUM))
	      yyalloc = YYSTACK_ALLOC_MAXIMUM;
	    if (yymsg != yymsgbuf)
	      YYSTACK_FREE (yymsg);
	    yymsg = (char *) YYSTACK_ALLOC (yyalloc);
	    if (yymsg)
	      yymsg_alloc = yyalloc;
	    else
	      {
		yymsg = yymsgbuf;
		yymsg_alloc = sizeof yymsgbuf;
	      }
	  }

	if (0 < yysize && yysize <= yymsg_alloc)
	  {
	    (void) yysyntax_error (yymsg, yystate, yychar);
	    yyerror (yymsg);
	  }
	else
	  {
	    yyerror (YY_("syntax error"));
	    if (yysize != 0)
	      goto yyexhaustedlab;
	  }
      }
#endif
    }



  if (yyerrstatus == 3)
    {
      /* If just tried and failed to reuse lookahead token after an
	 error, discard it.  */

      if (yychar <= YYEOF)
	{
	  /* Return failure if at end of input.  */
	  if (yychar == YYEOF)
	    YYABORT;
	}
      else
	{
	  yydestruct ("Error: discarding",
		      yytoken, &yylval);
	  yychar = YYEMPTY;
	}
    }

  /* Else will try to reuse lookahead token after shifting the error
     token.  */
  goto yyerrlab1;


/*---------------------------------------------------.
| yyerrorlab -- error raised explicitly by YYERROR.  |
`---------------------------------------------------*/
yyerrorlab:

  /* Pacify compilers like GCC when the user code never invokes
     YYERROR and the label yyerrorlab therefore never appears in user
     code.  */
  if (/*CONSTCOND*/ 0)
     goto yyerrorlab;

  /* Do not reclaim the symbols of the rule which action triggered
     this YYERROR.  */
  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);
  yystate = *yyssp;
  goto yyerrlab1;


/*-------------------------------------------------------------.
| yyerrlab1 -- common code for both syntax error and YYERROR.  |
`-------------------------------------------------------------*/
yyerrlab1:
  yyerrstatus = 3;	/* Each real token shifted decrements this.  */

  for (;;)
    {
      yyn = yypact[yystate];
      if (yyn != YYPACT_NINF)
	{
	  yyn += YYTERROR;
	  if (0 <= yyn && yyn <= YYLAST && yycheck[yyn] == YYTERROR)
	    {
	      yyn = yytable[yyn];
	      if (0 < yyn)
		break;
	    }
	}

      /* Pop the current state because it cannot handle the error token.  */
      if (yyssp == yyss)
	YYABORT;


      yydestruct ("Error: popping",
		  yystos[yystate], yyvsp);
      YYPOPSTACK (1);
      yystate = *yyssp;
      YY_STACK_PRINT (yyss, yyssp);
    }

  *++yyvsp = yylval;


  /* Shift the error token.  */
  YY_SYMBOL_PRINT ("Shifting", yystos[yyn], yyvsp, yylsp);

  yystate = yyn;
  goto yynewstate;


/*-------------------------------------.
| yyacceptlab -- YYACCEPT comes here.  |
`-------------------------------------*/
yyacceptlab:
  yyresult = 0;
  goto yyreturn;

/*-----------------------------------.
| yyabortlab -- YYABORT comes here.  |
`-----------------------------------*/
yyabortlab:
  yyresult = 1;
  goto yyreturn;

#if !defined(yyoverflow) || YYERROR_VERBOSE
/*-------------------------------------------------.
| yyexhaustedlab -- memory exhaustion comes here.  |
`-------------------------------------------------*/
yyexhaustedlab:
  yyerror (YY_("memory exhausted"));
  yyresult = 2;
  /* Fall through.  */
#endif

yyreturn:
  if (yychar != YYEMPTY)
     yydestruct ("Cleanup: discarding lookahead",
		 yytoken, &yylval);
  /* Do not reclaim the symbols of the rule which action triggered
     this YYABORT or YYACCEPT.  */
  YYPOPSTACK (yylen);
  YY_STACK_PRINT (yyss, yyssp);
  while (yyssp != yyss)
    {
      yydestruct ("Cleanup: popping",
		  yystos[*yyssp], yyvsp);
      YYPOPSTACK (1);
    }
#ifndef yyoverflow
  if (yyss != yyssa)
    YYSTACK_FREE (yyss);
#endif
#if YYERROR_VERBOSE
  if (yymsg != yymsgbuf)
    YYSTACK_FREE (yymsg);
#endif
  /* Make sure YYID is used.  */
  return YYID (yyresult);
}



/* Line 1675 of yacc.c  */
#line 623 "regparser.y"



regRegion* regParse( char* buf )
{
    // Needed to clear buff prior to parsing.
    regYYrestart(NULL);
  
    regRegion* regptr;
    char* ptr;

    // my_Gregion is declared externally
    my_Gregion = NULL;
    regParseStr = buf;
    
    // Needed to ensure extent is correctly set
    double fx[2] ={ -DBL_MAX, DBL_MAX };
    double fy[2] ={ -DBL_MAX, DBL_MAX };

    if ( !buf ) {
        return NULL;
    }
    
    ptr = buf;
    while( *ptr == ' ' || *ptr == '(' ) ptr++;
    if ( !isalpha( *ptr ) && *ptr != '!' )
    {
        // Not a region ( always begins with alpha or ! )
        return NULL;
    }

    regParseStrEnd = buf + strlen( buf );
    
    // uncomment to access debug mode
    // regYYdebug = 1;

    regYYparse();
    regptr = my_Gregion;

    // If we have successfully parsed a region then be sure to set the
    // appropriate bounds.
    if (regptr) {
        regExtent(regptr, fx, fy, regptr->xregbounds, regptr->yregbounds);
    }

    return regptr;
}


void regYYerror( char* msg )
{
    my_Gregion = NULL;
    return;
}


