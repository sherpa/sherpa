
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

/* All symbols defined below should begin with regYY or YY, to avoid
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



/* Copy the first part of user declarations.  */

/* Line 189 of yacc.c  */
#line 2 "regparser.y"

/*_C_INSERT_SAO_COPYRIGHT_HERE_(2007)_*/
/*_C_INSERT_GPL_LICENSE_HERE_*/
#include "region_priv.h"
#include "regpar.h"
int test_link_jcm( void );
#include <ctype.h>
extern char     *regParseStr;
extern char     *regParseStrEnd;
extern regRegion *my_Gregion;
static int world_coord;
static int world_size;


/* Line 189 of yacc.c  */
#line 88 "regparser.tab.c"

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
   enum regYYtokentype {
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

/* Line 214 of yacc.c  */
#line 23 "regparser.y"

  double dval;
  char str[1024];
  regRegion* my_region;
  regShape* my_shape;
  struct polyside
  { double *polyX; double *polyY; long polyS; } PolySide;



/* Line 214 of yacc.c  */
#line 181 "regparser.tab.c"
} YYSTYPE;
# define YYSTYPE_IS_TRIVIAL 1
# define regYYstype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
#endif


/* Copy the second part of user declarations.  */


/* Line 264 of yacc.c  */
#line 193 "regparser.tab.c"

#ifdef short
# undef short
#endif

#ifdef YYTYPE_UINT8
typedef YYTYPE_UINT8 regYYtype_uint8;
#else
typedef unsigned char regYYtype_uint8;
#endif

#ifdef YYTYPE_INT8
typedef YYTYPE_INT8 regYYtype_int8;
#elif (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
typedef signed char regYYtype_int8;
#else
typedef short int regYYtype_int8;
#endif

#ifdef YYTYPE_UINT16
typedef YYTYPE_UINT16 regYYtype_uint16;
#else
typedef unsigned short int regYYtype_uint16;
#endif

#ifdef YYTYPE_INT16
typedef YYTYPE_INT16 regYYtype_int16;
#else
typedef short int regYYtype_int16;
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
YYID (int regYYi)
#else
static int
YYID (regYYi)
    int regYYi;
#endif
{
  return regYYi;
}
#endif

#if ! defined regYYoverflow || YYERROR_VERBOSE

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
#endif /* ! defined regYYoverflow || YYERROR_VERBOSE */


#if (! defined regYYoverflow \
     && (! defined __cplusplus \
	 || (defined YYSTYPE_IS_TRIVIAL && YYSTYPE_IS_TRIVIAL)))

/* A type that is properly aligned for any stack member.  */
union regYYalloc
{
  regYYtype_int16 regYYss_alloc;
  YYSTYPE regYYvs_alloc;
};

/* The size of the maximum gap between one aligned stack and the next.  */
# define YYSTACK_GAP_MAXIMUM (sizeof (union regYYalloc) - 1)

/* The size of an array large to enough to hold all stacks, each with
   N elements.  */
# define YYSTACK_BYTES(N) \
     ((N) * (sizeof (regYYtype_int16) + sizeof (YYSTYPE)) \
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
	  YYSIZE_T regYYi;				\
	  for (regYYi = 0; regYYi < (Count); regYYi++)	\
	    (To)[regYYi] = (From)[regYYi];		\
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
	YYSIZE_T regYYnewbytes;						\
	YYCOPY (&regYYptr->Stack_alloc, Stack, regYYsize);			\
	Stack = &regYYptr->Stack_alloc;					\
	regYYnewbytes = regYYstacksize * sizeof (*Stack) + YYSTACK_GAP_MAXIMUM; \
	regYYptr += regYYnewbytes / sizeof (*regYYptr);				\
      }									\
    while (YYID (0))

#endif

/* YYFINAL -- State number of the termination state.  */
#define YYFINAL  46
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   375

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  34
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  9
/* YYNRULES -- Number of rules.  */
#define YYNRULES  57
/* YYNRULES -- Number of states.  */
#define YYNSTATES  289

/* YYTRANSLATE(YYLEX) -- Bison symbol number corresponding to YYLEX.  */
#define YYUNDEFTOK  2
#define YYMAXUTOK   279

#define YYTRANSLATE(YYX)						\
  ((unsigned int) (YYX) <= YYMAXUTOK ? regYYtranslate[YYX] : YYUNDEFTOK)

/* YYTRANSLATE[YYLEX] -- Bison symbol number corresponding to YYLEX.  */
static const regYYtype_uint8 regYYtranslate[] =
{
       0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,    28,     2,     2,     2,     2,    27,
      31,    32,     2,     2,    33,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,    26,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
      25,     2,     2,     2,     2,    30,     2,     2,     2,     2,
       2,     2,    29,     2,     2,     2,     2,     2,     2,     2,
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
      15,    16,    17,    18,    19,    20,    21,    22,    23,    24
};

#if YYDEBUG
/* YYPRHS[YYN] -- Index of the first RHS symbol of rule number YYN in
   YYRHS.  */
static const regYYtype_uint16 regYYprhs[] =
{
       0,     0,     3,     7,    10,    14,    18,    21,    23,    26,
      28,    31,    37,    39,    42,    48,    55,    57,    60,    62,
      65,    68,    71,    73,    76,    79,    82,    85,    89,    92,
      96,   101,   110,   119,   129,   140,   152,   163,   175,   188,
     202,   215,   229,   244,   260,   271,   283,   290,   298,   309,
     321,   332,   344,   357,   371,   376,   382,   388
};

/* YYRHS -- A `-1'-separated list of the rules' RHS.  */
static const regYYtype_int8 regYYrhs[] =
{
      35,     0,    -1,    35,    16,    41,    -1,    35,    41,    -1,
      35,    15,    41,    -1,    35,    17,    41,    -1,    16,    41,
      -1,    41,    -1,     1,    24,    -1,    38,    -1,    38,    25,
      -1,    22,    26,    22,    26,    22,    -1,    38,    -1,    38,
      25,    -1,    22,    26,    22,    26,    22,    -1,    17,    22,
      26,    22,    26,    22,    -1,    22,    -1,    17,    22,    -1,
      38,    -1,    38,    27,    -1,    38,    28,    -1,    38,    25,
      -1,    38,    -1,    38,    29,    -1,    38,    30,    -1,    38,
      25,    -1,    38,    27,    -1,    38,    27,    27,    -1,    38,
      28,    -1,    18,    31,    32,    -1,    14,    18,    31,    32,
      -1,    19,    31,    36,    33,    37,    33,    20,    32,    -1,
       3,    31,    36,    33,    37,    33,    40,    32,    -1,    14,
       3,    31,    36,    33,    37,    33,    40,    32,    -1,    10,
      31,    36,    33,    37,    33,    40,    33,    40,    32,    -1,
      14,    10,    31,    36,    33,    37,    33,    40,    33,    40,
      32,    -1,     4,    31,    36,    33,    37,    33,    40,    33,
      40,    32,    -1,    14,     4,    31,    36,    33,    37,    33,
      40,    33,    40,    32,    -1,     4,    31,    36,    33,    37,
      33,    40,    33,    40,    33,    39,    32,    -1,    14,     4,
      31,    36,    33,    37,    33,    40,    33,    40,    33,    39,
      32,    -1,     9,    31,    36,    33,    37,    33,    40,    33,
      40,    33,    39,    32,    -1,    14,     9,    31,    36,    33,
      37,    33,    40,    33,    40,    33,    39,    32,    -1,     6,
      31,    36,    33,    37,    33,    40,    33,    40,    33,    39,
      33,    39,    32,    -1,    14,     6,    31,    36,    33,    37,
      33,    40,    33,    40,    33,    39,    33,    39,    32,    -1,
      11,    31,    36,    33,    37,    33,    39,    33,    39,    32,
      -1,    14,    11,    31,    36,    33,    37,    33,    39,    33,
      39,    32,    -1,    12,    31,    36,    33,    37,    32,    -1,
      14,    12,    31,    36,    33,    37,    32,    -1,     7,    31,
      36,    33,    37,    33,    36,    33,    37,    32,    -1,    14,
       7,    31,    36,    33,    37,    33,    36,    33,    37,    32,
      -1,    21,    31,    36,    33,    37,    33,    36,    33,    37,
      32,    -1,    14,    21,    31,    36,    33,    37,    33,    36,
      33,    37,    32,    -1,     8,    31,    36,    33,    37,    33,
      40,    33,    40,    33,    39,    32,    -1,    14,     8,    31,
      36,    33,    37,    33,    40,    33,    40,    33,    39,    32,
      -1,     5,    31,    42,    32,    -1,    14,     5,    31,    42,
      32,    -1,    42,    33,    36,    33,    37,    -1,    36,    33,
      37,    -1
};

/* YYRLINE[YYN] -- source line where rule number YYN was defined.  */
static const regYYtype_uint16 regYYrline[] =
{
       0,    45,    45,    46,    47,    48,    49,    54,    59,    65,
      67,    69,    80,    82,    84,    92,   103,   105,   110,   112,
     114,   116,   121,   123,   125,   127,   129,   131,   133,   139,
     147,   155,   166,   182,   197,   218,   239,   255,   270,   285,
     300,   315,   330,   349,   368,   378,   388,   398,   408,   418,
     428,   439,   449,   464,   479,   491,   506,   515
};
#endif

#if YYDEBUG || YYERROR_VERBOSE || YYTOKEN_TABLE
/* YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals.  */
static const char *const regYYtname[] =
{
  "$end", "error", "$undefined", "REG_CIRCLE", "REG_BOX", "REG_POLY",
  "REG_PIE", "REG_RECT", "REG_ELL", "REG_ROTBOX", "REG_ANNULUS",
  "REG_SECTOR", "REG_POINT", "REG_ROTRECT", "REG_NEG", "REG_AND", "REG_OR",
  "REG_MINUS", "REG_FIELD", "REG_TEXT", "REG_STRING", "REG_LINE",
  "REG_NUMBER", "colName", "\"\"", "'d'", "':'", "'\\''", "'\"'", "'p'",
  "'i'", "'('", "')'", "','", "$accept", "reg_comp", "reg_xcoord",
  "reg_ycoord", "reg_value", "reg_angle", "reg_size", "reg_shape",
  "reg_ordered_pair", 0
};
#endif

# ifdef YYPRINT
/* YYTOKNUM[YYLEX-NUM] -- Internal token number corresponding to
   token YYLEX-NUM.  */
static const regYYtype_uint16 regYYtoknum[] =
{
       0,   256,   257,   258,   259,   260,   261,   262,   263,   264,
     265,   266,   267,   268,   269,   270,   271,   272,   273,   274,
     275,   276,   277,   278,   279,   100,    58,    39,    34,   112,
     105,    40,    41,    44
};
# endif

/* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const regYYtype_uint8 regYYr1[] =
{
       0,    34,    35,    35,    35,    35,    35,    35,    35,    36,
      36,    36,    37,    37,    37,    37,    38,    38,    39,    39,
      39,    39,    40,    40,    40,    40,    40,    40,    40,    41,
      41,    41,    41,    41,    41,    41,    41,    41,    41,    41,
      41,    41,    41,    41,    41,    41,    41,    41,    41,    41,
      41,    41,    41,    41,    41,    41,    42,    42
};

/* YYR2[YYN] -- Number of symbols composing right hand side of rule YYN.  */
static const regYYtype_uint8 regYYr2[] =
{
       0,     2,     3,     2,     3,     3,     2,     1,     2,     1,
       2,     5,     1,     2,     5,     6,     1,     2,     1,     2,
       2,     2,     1,     2,     2,     2,     2,     3,     2,     3,
       4,     8,     8,     9,    10,    11,    10,    11,    12,    13,
      12,    13,    14,    15,    10,    11,     6,     7,    10,    11,
      10,    11,    12,    13,     4,     5,     5,     3
};

/* YYDEFACT[STATE-NAME] -- Default rule to reduce with in state
   STATE-NUM when YYTABLE doesn't specify something else to do.  Zero
   means the default is an error.  */
static const regYYtype_uint8 regYYdefact[] =
{
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     7,     8,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     6,     0,     0,     0,     1,     0,     0,     0,
       3,     0,    16,     0,     9,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,    29,     0,     0,
       4,     2,     5,    17,     0,     0,    10,     0,     0,    54,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,    30,     0,
       0,     0,     0,     0,    16,     0,    12,     0,    57,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,    55,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,    17,     0,     0,    13,     0,     0,     0,     0,     0,
       0,     0,     0,    46,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,    11,     0,     0,    16,
      22,     0,     0,    56,     0,     0,     0,     0,     0,    18,
       0,     0,     0,     0,     0,     0,     0,     0,     0,    47,
       0,     0,     0,     0,     0,    25,    26,    28,    23,    24,
      32,     0,     0,     0,     0,     0,     0,    21,    19,    20,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
      31,     0,     0,    14,    27,     0,     0,     0,     0,     0,
       0,     0,    33,     0,     0,     0,     0,     0,     0,     0,
       0,     0,    15,    36,     0,     0,    48,     0,     0,    34,
      44,     0,     0,     0,     0,     0,     0,     0,     0,    50,
       0,     0,     0,     0,    37,     0,     0,    49,     0,     0,
      35,    45,    51,    38,     0,    52,    40,     0,     0,     0,
       0,     0,    39,     0,    53,    41,    42,     0,    43
};

/* YYDEFGOTO[NTERM-NUM].  */
static const regYYtype_int16 regYYdefgoto[] =
{
      -1,    17,    56,   115,    54,   180,   171,    18,    57
};

/* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
   STATE-NUM.  */
#define YYPACT_NINF -102
static const regYYtype_int16 regYYpact[] =
{
     310,   -23,    -4,    32,    38,    45,    53,    55,    73,    75,
      76,    77,    10,   329,    80,    81,    82,   273,  -102,  -102,
      22,    22,    22,    22,    22,    22,    22,    22,    22,    22,
      83,    84,    91,    92,    93,    95,    96,    97,    98,    99,
     100,   101,  -102,   102,    22,    22,  -102,   329,   329,   329,
    -102,   111,   113,   104,   110,   107,   108,   -28,   112,   123,
     125,   128,   133,   136,   137,    22,    22,    22,    22,    22,
      22,    22,    22,    22,    22,   139,    22,  -102,   142,   143,
    -102,  -102,  -102,  -102,   120,    71,  -102,    71,    71,  -102,
      22,    71,    71,    71,    71,    71,    71,    71,   144,   145,
      -9,   146,   152,   153,   154,   155,   157,   158,  -102,   159,
      71,    71,   168,   150,   169,   174,   183,   176,  -102,   177,
     179,   181,   185,   188,   189,   190,   192,    71,    71,  -102,
      71,    71,    71,    71,    71,    71,    71,    71,   193,   194,
     203,   202,   207,    88,  -102,    88,    71,    88,    22,    88,
      88,    88,    88,  -102,   197,   198,   199,   200,   201,   205,
     217,   221,   224,   225,   215,    22,  -102,   235,   233,  -102,
      43,   228,   229,  -102,   230,   231,   232,   236,   241,    67,
     242,    88,    88,    88,    22,    88,    88,    88,    88,  -102,
      22,   234,   253,   267,   239,  -102,   268,  -102,  -102,  -102,
    -102,    88,    88,    71,    88,    88,    88,  -102,  -102,  -102,
      88,   264,   265,   272,   274,   275,   276,   277,   290,   292,
    -102,    71,   279,  -102,  -102,     8,   294,   313,   309,   318,
     314,   320,  -102,    88,    88,    71,    88,    88,    88,    88,
      71,   321,  -102,  -102,    88,    88,  -102,    88,    88,  -102,
    -102,    24,   322,   324,   325,   326,   328,   330,   331,  -102,
     332,   333,   335,   336,  -102,    88,    88,  -102,    88,    88,
    -102,  -102,  -102,  -102,    88,  -102,  -102,   337,   338,   340,
     341,   342,  -102,    88,  -102,  -102,  -102,   343,  -102
};

/* YYPGOTO[NTERM-NUM].  */
static const regYYtype_int16 regYYpgoto[] =
{
    -102,  -102,     9,   109,   -85,  -101,    66,    42,   287
};

/* YYTABLE[YYPACT[STATE-NUM]].  What to do in state STATE-NUM.  If
   positive, shift that token.  If negative, reduce the rule which
   number is the opposite.  If zero, do what YYDEFACT says.
   If YYTABLE_NINF, syntax error.  */
#define YYTABLE_NINF -1
static const regYYtype_uint16 regYYtable[] =
{
     116,    19,   116,   116,    89,    90,   116,   116,   116,   116,
     116,   116,   116,    30,    31,    32,    33,    34,    35,    36,
      37,    38,    39,   129,    90,   116,   116,    20,    40,    53,
      55,    41,    58,    59,    60,    61,    62,    63,    64,    51,
     243,   244,   116,   116,    52,   116,   116,   116,   116,   116,
     116,   116,   116,    78,    79,    42,   264,   265,   170,    50,
     170,   116,   170,    21,   170,   170,   170,   179,   195,    22,
     196,   197,   198,   199,    98,    99,    23,   101,   102,   103,
     104,   105,   106,   107,    24,   109,    25,   218,   113,    80,
      81,    82,   207,   114,   208,   209,   170,   170,   170,   119,
     170,   170,   170,   179,    26,    51,    27,    28,    29,   231,
     169,    43,    44,    45,    65,    66,   170,   170,   116,   170,
     170,   170,    67,    68,    69,   179,    70,    71,    72,    73,
      74,    75,    76,    83,    77,    86,   116,    85,   257,    84,
      87,    88,   112,   260,   261,    91,   262,   263,   170,   170,
     116,   170,   170,   170,   179,   116,    92,   175,    93,   179,
     179,    94,   179,   179,   277,   278,    95,   279,   280,    96,
      97,   108,   141,   281,   192,   110,   111,   127,   128,   130,
     179,   179,   287,   179,   179,   131,   132,   133,   134,   179,
     135,   136,   137,   214,   140,   142,   117,   118,   179,   219,
     120,   121,   122,   123,   124,   125,   126,   143,   144,   145,
     146,   172,   147,   174,   148,   176,   177,   178,   149,   138,
     139,   150,   151,   152,   153,   166,   164,   165,   167,   168,
     181,   182,   183,   184,   185,   191,   154,   155,   186,   156,
     157,   158,   159,   160,   161,   162,   163,   211,   212,   213,
     187,   215,   216,   217,   188,   173,   189,   193,   190,   194,
     200,   223,   201,   202,   203,   204,   220,   225,   226,   205,
     228,   229,   230,    46,   206,   210,     2,     3,     4,     5,
       6,     7,     8,     9,    10,    11,   221,    12,    47,    48,
      49,    14,    15,   222,    16,   224,   232,     0,   233,   251,
     252,   242,   254,   255,   256,   234,     0,   235,   236,   237,
     238,     1,   227,     2,     3,     4,     5,     6,     7,     8,
       9,    10,    11,   239,    12,   240,    13,   245,    14,    15,
     241,    16,     2,     3,     4,     5,     6,     7,     8,     9,
      10,    11,   247,    12,   253,   246,   249,    14,    15,   258,
      16,   248,   250,   259,   100,   266,   267,     0,   268,   269,
     270,     0,   271,   272,   273,     0,   274,   275,   276,   282,
       0,   283,   284,   285,   286,   288
};

static const regYYtype_int16 regYYcheck[] =
{
      85,    24,    87,    88,    32,    33,    91,    92,    93,    94,
      95,    96,    97,     3,     4,     5,     6,     7,     8,     9,
      10,    11,    12,    32,    33,   110,   111,    31,    18,    20,
      21,    21,    23,    24,    25,    26,    27,    28,    29,    17,
      32,    33,   127,   128,    22,   130,   131,   132,   133,   134,
     135,   136,   137,    44,    45,    13,    32,    33,   143,    17,
     145,   146,   147,    31,   149,   150,   151,   152,    25,    31,
      27,    28,    29,    30,    65,    66,    31,    68,    69,    70,
      71,    72,    73,    74,    31,    76,    31,   188,    17,    47,
      48,    49,    25,    22,    27,    28,   181,   182,   183,    90,
     185,   186,   187,   188,    31,    17,    31,    31,    31,   210,
      22,    31,    31,    31,    31,    31,   201,   202,   203,   204,
     205,   206,    31,    31,    31,   210,    31,    31,    31,    31,
      31,    31,    31,    22,    32,    25,   221,    33,   239,    26,
      33,    33,    22,   244,   245,    33,   247,   248,   233,   234,
     235,   236,   237,   238,   239,   240,    33,   148,    33,   244,
     245,    33,   247,   248,   265,   266,    33,   268,   269,    33,
      33,    32,    22,   274,   165,    33,    33,    33,    33,    33,
     265,   266,   283,   268,   269,    33,    33,    33,    33,   274,
      33,    33,    33,   184,    26,    26,    87,    88,   283,   190,
      91,    92,    93,    94,    95,    96,    97,    33,    25,    33,
      33,   145,    33,   147,    33,   149,   150,   151,    33,   110,
     111,    33,    33,    33,    32,    22,    33,    33,    26,    22,
      33,    33,    33,    33,    33,    20,   127,   128,    33,   130,
     131,   132,   133,   134,   135,   136,   137,   181,   182,   183,
      33,   185,   186,   187,    33,   146,    32,    22,    33,    26,
      32,    22,    33,    33,    33,    33,    32,   201,   202,    33,
     204,   205,   206,     0,    33,    33,     3,     4,     5,     6,
       7,     8,     9,    10,    11,    12,    33,    14,    15,    16,
      17,    18,    19,    26,    21,    27,    32,    -1,    33,   233,
     234,    22,   236,   237,   238,    33,    -1,    33,    33,    33,
      33,     1,   203,     3,     4,     5,     6,     7,     8,     9,
      10,    11,    12,    33,    14,    33,    16,    33,    18,    19,
     221,    21,     3,     4,     5,     6,     7,     8,     9,    10,
      11,    12,    33,    14,   235,    32,    32,    18,    19,   240,
      21,    33,    32,    32,    67,    33,    32,    -1,    33,    33,
      32,    -1,    32,    32,    32,    -1,    33,    32,    32,    32,
      -1,    33,    32,    32,    32,    32
};

/* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
   symbol of state STATE-NUM.  */
static const regYYtype_uint8 regYYstos[] =
{
       0,     1,     3,     4,     5,     6,     7,     8,     9,    10,
      11,    12,    14,    16,    18,    19,    21,    35,    41,    24,
      31,    31,    31,    31,    31,    31,    31,    31,    31,    31,
       3,     4,     5,     6,     7,     8,     9,    10,    11,    12,
      18,    21,    41,    31,    31,    31,     0,    15,    16,    17,
      41,    17,    22,    36,    38,    36,    36,    42,    36,    36,
      36,    36,    36,    36,    36,    31,    31,    31,    31,    31,
      31,    31,    31,    31,    31,    31,    31,    32,    36,    36,
      41,    41,    41,    22,    26,    33,    25,    33,    33,    32,
      33,    33,    33,    33,    33,    33,    33,    33,    36,    36,
      42,    36,    36,    36,    36,    36,    36,    36,    32,    36,
      33,    33,    22,    17,    22,    37,    38,    37,    37,    36,
      37,    37,    37,    37,    37,    37,    37,    33,    33,    32,
      33,    33,    33,    33,    33,    33,    33,    33,    37,    37,
      26,    22,    26,    33,    25,    33,    33,    33,    33,    33,
      33,    33,    33,    32,    37,    37,    37,    37,    37,    37,
      37,    37,    37,    37,    33,    33,    22,    26,    22,    22,
      38,    40,    40,    37,    40,    36,    40,    40,    40,    38,
      39,    33,    33,    33,    33,    33,    33,    33,    33,    32,
      33,    20,    36,    22,    26,    25,    27,    28,    29,    30,
      32,    33,    33,    33,    33,    33,    33,    25,    27,    28,
      33,    40,    40,    40,    36,    40,    40,    40,    39,    36,
      32,    33,    26,    22,    27,    40,    40,    37,    40,    40,
      40,    39,    32,    33,    33,    33,    33,    33,    33,    33,
      33,    37,    22,    32,    33,    33,    32,    33,    33,    32,
      32,    40,    40,    37,    40,    40,    40,    39,    37,    32,
      39,    39,    39,    39,    32,    33,    33,    32,    33,    33,
      32,    32,    32,    32,    33,    32,    32,    39,    39,    39,
      39,    39,    32,    33,    32,    32,    32,    39,    32
};

#define regYYerrok		(regYYerrstatus = 0)
#define regYYclearin	(regYYchar = YYEMPTY)
#define YYEMPTY		(-2)
#define YYEOF		0

#define YYACCEPT	goto regYYacceptlab
#define YYABORT		goto regYYabortlab
#define YYERROR		goto regYYerrorlab


/* Like YYERROR except do call regYYerror.  This remains here temporarily
   to ease the transition to the new meaning of YYERROR, for GCC.
   Once GCC version 2 has supplanted version 1, this can go.  */

#define YYFAIL		goto regYYerrlab

#define YYRECOVERING()  (!!regYYerrstatus)

#define YYBACKUP(Token, Value)					\
do								\
  if (regYYchar == YYEMPTY && regYYlen == 1)				\
    {								\
      regYYchar = (Token);						\
      regYYlval = (Value);						\
      regYYtoken = YYTRANSLATE (regYYchar);				\
      YYPOPSTACK (1);						\
      goto regYYbackup;						\
    }								\
  else								\
    {								\
      regYYerror (YY_("syntax error: cannot back up")); \
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


/* YYLEX -- calling `regYYlex' with the right arguments.  */

#ifdef YYLEX_PARAM
# define YYLEX regYYlex (YYLEX_PARAM)
#else
# define YYLEX regYYlex ()
#endif

/* Enable debugging if requested.  */
#if YYDEBUG

# ifndef YYFPRINTF
#  include <stdio.h> /* INFRINGES ON USER NAME SPACE */
#  define YYFPRINTF fprintf
# endif

# define YYDPRINTF(Args)			\
do {						\
  if (regYYdebug)					\
    YYFPRINTF Args;				\
} while (YYID (0))

# define YY_SYMBOL_PRINT(Title, Type, Value, Location)			  \
do {									  \
  if (regYYdebug)								  \
    {									  \
      YYFPRINTF (stderr, "%s ", Title);					  \
      regYY_symbol_print (stderr,						  \
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
regYY_symbol_value_print (FILE *regYYoutput, int regYYtype, YYSTYPE const * const regYYvaluep)
#else
static void
regYY_symbol_value_print (regYYoutput, regYYtype, regYYvaluep)
    FILE *regYYoutput;
    int regYYtype;
    YYSTYPE const * const regYYvaluep;
#endif
{
  if (!regYYvaluep)
    return;
# ifdef YYPRINT
  if (regYYtype < YYNTOKENS)
    YYPRINT (regYYoutput, regYYtoknum[regYYtype], *regYYvaluep);
# else
  YYUSE (regYYoutput);
# endif
  switch (regYYtype)
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
regYY_symbol_print (FILE *regYYoutput, int regYYtype, YYSTYPE const * const regYYvaluep)
#else
static void
regYY_symbol_print (regYYoutput, regYYtype, regYYvaluep)
    FILE *regYYoutput;
    int regYYtype;
    YYSTYPE const * const regYYvaluep;
#endif
{
  if (regYYtype < YYNTOKENS)
    YYFPRINTF (regYYoutput, "token %s (", regYYtname[regYYtype]);
  else
    YYFPRINTF (regYYoutput, "nterm %s (", regYYtname[regYYtype]);

  regYY_symbol_value_print (regYYoutput, regYYtype, regYYvaluep);
  YYFPRINTF (regYYoutput, ")");
}

/*------------------------------------------------------------------.
| regYY_stack_print -- Print the state stack from its BOTTOM up to its |
| TOP (included).                                                   |
`------------------------------------------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
regYY_stack_print (regYYtype_int16 *regYYbottom, regYYtype_int16 *regYYtop)
#else
static void
regYY_stack_print (regYYbottom, regYYtop)
    regYYtype_int16 *regYYbottom;
    regYYtype_int16 *regYYtop;
#endif
{
  YYFPRINTF (stderr, "Stack now");
  for (; regYYbottom <= regYYtop; regYYbottom++)
    {
      int regYYbot = *regYYbottom;
      YYFPRINTF (stderr, " %d", regYYbot);
    }
  YYFPRINTF (stderr, "\n");
}

# define YY_STACK_PRINT(Bottom, Top)				\
do {								\
  if (regYYdebug)							\
    regYY_stack_print ((Bottom), (Top));				\
} while (YYID (0))


/*------------------------------------------------.
| Report that the YYRULE is going to be reduced.  |
`------------------------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
regYY_reduce_print (YYSTYPE *regYYvsp, int regYYrule)
#else
static void
regYY_reduce_print (regYYvsp, regYYrule)
    YYSTYPE *regYYvsp;
    int regYYrule;
#endif
{
  int regYYnrhs = regYYr2[regYYrule];
  int regYYi;
  unsigned long int regYYlno = regYYrline[regYYrule];
  YYFPRINTF (stderr, "Reducing stack by rule %d (line %lu):\n",
	     regYYrule - 1, regYYlno);
  /* The symbols being reduced.  */
  for (regYYi = 0; regYYi < regYYnrhs; regYYi++)
    {
      YYFPRINTF (stderr, "   $%d = ", regYYi + 1);
      regYY_symbol_print (stderr, regYYrhs[regYYprhs[regYYrule] + regYYi],
		       &(regYYvsp[(regYYi + 1) - (regYYnrhs)])
		       		       );
      YYFPRINTF (stderr, "\n");
    }
}

# define YY_REDUCE_PRINT(Rule)		\
do {					\
  if (regYYdebug)				\
    regYY_reduce_print (regYYvsp, Rule); \
} while (YYID (0))

/* Nonzero means print parse trace.  It is left uninitialized so that
   multiple parsers can coexist.  */
int regYYdebug;
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

# ifndef regYYstrlen
#  if defined __GLIBC__ && defined _STRING_H
#   define regYYstrlen strlen
#  else
/* Return the length of YYSTR.  */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static YYSIZE_T
regYYstrlen (const char *regYYstr)
#else
static YYSIZE_T
regYYstrlen (regYYstr)
    const char *regYYstr;
#endif
{
  YYSIZE_T regYYlen;
  for (regYYlen = 0; regYYstr[regYYlen]; regYYlen++)
    continue;
  return regYYlen;
}
#  endif
# endif

# ifndef regYYstpcpy
#  if defined __GLIBC__ && defined _STRING_H && defined _GNU_SOURCE
#   define regYYstpcpy stpcpy
#  else
/* Copy YYSRC to YYDEST, returning the address of the terminating '\0' in
   YYDEST.  */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static char *
regYYstpcpy (char *regYYdest, const char *regYYsrc)
#else
static char *
regYYstpcpy (regYYdest, regYYsrc)
    char *regYYdest;
    const char *regYYsrc;
#endif
{
  char *regYYd = regYYdest;
  const char *regYYs = regYYsrc;

  while ((*regYYd++ = *regYYs++) != '\0')
    continue;

  return regYYd - 1;
}
#  endif
# endif

# ifndef regYYtnamerr
/* Copy to YYRES the contents of YYSTR after stripping away unnecessary
   quotes and backslashes, so that it's suitable for regYYerror.  The
   heuristic is that double-quoting is unnecessary unless the string
   contains an apostrophe, a comma, or backslash (other than
   backslash-backslash).  YYSTR is taken from regYYtname.  If YYRES is
   null, do not copy; instead, return the length of what the result
   would have been.  */
static YYSIZE_T
regYYtnamerr (char *regYYres, const char *regYYstr)
{
  if (*regYYstr == '"')
    {
      YYSIZE_T regYYn = 0;
      char const *regYYp = regYYstr;

      for (;;)
	switch (*++regYYp)
	  {
	  case '\'':
	  case ',':
	    goto do_not_strip_quotes;

	  case '\\':
	    if (*++regYYp != '\\')
	      goto do_not_strip_quotes;
	    /* Fall through.  */
	  default:
	    if (regYYres)
	      regYYres[regYYn] = *regYYp;
	    regYYn++;
	    break;

	  case '"':
	    if (regYYres)
	      regYYres[regYYn] = '\0';
	    return regYYn;
	  }
    do_not_strip_quotes: ;
    }

  if (! regYYres)
    return regYYstrlen (regYYstr);

  return regYYstpcpy (regYYres, regYYstr) - regYYres;
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
regYYsyntax_error (char *regYYresult, int regYYstate, int regYYchar)
{
  int regYYn = regYYpact[regYYstate];

  if (! (YYPACT_NINF < regYYn && regYYn <= YYLAST))
    return 0;
  else
    {
      int regYYtype = YYTRANSLATE (regYYchar);
      YYSIZE_T regYYsize0 = regYYtnamerr (0, regYYtname[regYYtype]);
      YYSIZE_T regYYsize = regYYsize0;
      YYSIZE_T regYYsize1;
      int regYYsize_overflow = 0;
      enum { YYERROR_VERBOSE_ARGS_MAXIMUM = 5 };
      char const *regYYarg[YYERROR_VERBOSE_ARGS_MAXIMUM];
      int regYYx;

# if 0
      /* This is so xgettext sees the translatable formats that are
	 constructed on the fly.  */
      YY_("syntax error, unexpected %s");
      YY_("syntax error, unexpected %s, expecting %s");
      YY_("syntax error, unexpected %s, expecting %s or %s");
      YY_("syntax error, unexpected %s, expecting %s or %s or %s");
      YY_("syntax error, unexpected %s, expecting %s or %s or %s or %s");
# endif
      char *regYYfmt;
      char const *regYYf;
      static char const regYYunexpected[] = "syntax error, unexpected %s";
      static char const regYYexpecting[] = ", expecting %s";
      static char const regYYor[] = " or %s";
      char regYYformat[sizeof regYYunexpected
		    + sizeof regYYexpecting - 1
		    + ((YYERROR_VERBOSE_ARGS_MAXIMUM - 2)
		       * (sizeof regYYor - 1))];
      char const *regYYprefix = regYYexpecting;

      /* Start YYX at -YYN if negative to avoid negative indexes in
	 YYCHECK.  */
      int regYYxbegin = regYYn < 0 ? -regYYn : 0;

      /* Stay within bounds of both regYYcheck and regYYtname.  */
      int regYYchecklim = YYLAST - regYYn + 1;
      int regYYxend = regYYchecklim < YYNTOKENS ? regYYchecklim : YYNTOKENS;
      int regYYcount = 1;

      regYYarg[0] = regYYtname[regYYtype];
      regYYfmt = regYYstpcpy (regYYformat, regYYunexpected);

      for (regYYx = regYYxbegin; regYYx < regYYxend; ++regYYx)
	if (regYYcheck[regYYx + regYYn] == regYYx && regYYx != YYTERROR)
	  {
	    if (regYYcount == YYERROR_VERBOSE_ARGS_MAXIMUM)
	      {
		regYYcount = 1;
		regYYsize = regYYsize0;
		regYYformat[sizeof regYYunexpected - 1] = '\0';
		break;
	      }
	    regYYarg[regYYcount++] = regYYtname[regYYx];
	    regYYsize1 = regYYsize + regYYtnamerr (0, regYYtname[regYYx]);
	    regYYsize_overflow |= (regYYsize1 < regYYsize);
	    regYYsize = regYYsize1;
	    regYYfmt = regYYstpcpy (regYYfmt, regYYprefix);
	    regYYprefix = regYYor;
	  }

      regYYf = YY_(regYYformat);
      regYYsize1 = regYYsize + regYYstrlen (regYYf);
      regYYsize_overflow |= (regYYsize1 < regYYsize);
      regYYsize = regYYsize1;

      if (regYYsize_overflow)
	return YYSIZE_MAXIMUM;

      if (regYYresult)
	{
	  /* Avoid sprintf, as that infringes on the user's name space.
	     Don't have undefined behavior even if the translation
	     produced a string with the wrong number of "%s"s.  */
	  char *regYYp = regYYresult;
	  int regYYi = 0;
	  while ((*regYYp = *regYYf) != '\0')
	    {
	      if (*regYYp == '%' && regYYf[1] == 's' && regYYi < regYYcount)
		{
		  regYYp += regYYtnamerr (regYYp, regYYarg[regYYi++]);
		  regYYf += 2;
		}
	      else
		{
		  regYYp++;
		  regYYf++;
		}
	    }
	}
      return regYYsize;
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
regYYdestruct (const char *regYYmsg, int regYYtype, YYSTYPE *regYYvaluep)
#else
static void
regYYdestruct (regYYmsg, regYYtype, regYYvaluep)
    const char *regYYmsg;
    int regYYtype;
    YYSTYPE *regYYvaluep;
#endif
{
  YYUSE (regYYvaluep);

  if (!regYYmsg)
    regYYmsg = "Deleting";
  YY_SYMBOL_PRINT (regYYmsg, regYYtype, regYYvaluep, regYYlocationp);

  switch (regYYtype)
    {

      default:
	break;
    }
}

/* Prevent warnings from -Wmissing-prototypes.  */
#ifdef YYPARSE_PARAM
#if defined __STDC__ || defined __cplusplus
int regYYparse (void *YYPARSE_PARAM);
#else
int regYYparse ();
#endif
#else /* ! YYPARSE_PARAM */
#if defined __STDC__ || defined __cplusplus
int regYYparse (void);
#else
int regYYparse ();
#endif
#endif /* ! YYPARSE_PARAM */


/* The lookahead symbol.  */
int regYYchar;

/* The semantic value of the lookahead symbol.  */
YYSTYPE regYYlval;

/* Number of syntax errors so far.  */
int regYYnerrs;



/*-------------------------.
| regYYparse or regYYpush_parse.  |
`-------------------------*/

#ifdef YYPARSE_PARAM
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
int
regYYparse (void *YYPARSE_PARAM)
#else
int
regYYparse (YYPARSE_PARAM)
    void *YYPARSE_PARAM;
#endif
#else /* ! YYPARSE_PARAM */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
int
regYYparse (void)
#else
int
regYYparse ()

#endif
#endif
{


    int regYYstate;
    /* Number of tokens to shift before error messages enabled.  */
    int regYYerrstatus;

    /* The stacks and their tools:
       `regYYss': related to states.
       `regYYvs': related to semantic values.

       Refer to the stacks thru separate pointers, to allow regYYoverflow
       to reallocate them elsewhere.  */

    /* The state stack.  */
    regYYtype_int16 regYYssa[YYINITDEPTH];
    regYYtype_int16 *regYYss;
    regYYtype_int16 *regYYssp;

    /* The semantic value stack.  */
    YYSTYPE regYYvsa[YYINITDEPTH];
    YYSTYPE *regYYvs;
    YYSTYPE *regYYvsp;

    YYSIZE_T regYYstacksize;

  int regYYn;
  int regYYresult;
  /* Lookahead token as an internal (translated) token number.  */
  int regYYtoken;
  /* The variables used to return semantic value and location from the
     action routines.  */
  YYSTYPE regYYval;

#if YYERROR_VERBOSE
  /* Buffer for error messages, and its allocated size.  */
  char regYYmsgbuf[128];
  char *regYYmsg = regYYmsgbuf;
  YYSIZE_T regYYmsg_alloc = sizeof regYYmsgbuf;
#endif

#define YYPOPSTACK(N)   (regYYvsp -= (N), regYYssp -= (N))

  /* The number of symbols on the RHS of the reduced rule.
     Keep to zero when no symbol should be popped.  */
  int regYYlen = 0;

  regYYtoken = 0;
  regYYss = regYYssa;
  regYYvs = regYYvsa;
  regYYstacksize = YYINITDEPTH;

  YYDPRINTF ((stderr, "Starting parse\n"));

  regYYstate = 0;
  regYYerrstatus = 0;
  regYYnerrs = 0;
  regYYchar = YYEMPTY; /* Cause a token to be read.  */

  /* Initialize stack pointers.
     Waste one element of value and location stack
     so that they stay on the same level as the state stack.
     The wasted elements are never initialized.  */
  regYYssp = regYYss;
  regYYvsp = regYYvs;

  goto regYYsetstate;

/*------------------------------------------------------------.
| regYYnewstate -- Push a new state, which is found in regYYstate.  |
`------------------------------------------------------------*/
 regYYnewstate:
  /* In all cases, when you get here, the value and location stacks
     have just been pushed.  So pushing a state here evens the stacks.  */
  regYYssp++;

 regYYsetstate:
  *regYYssp = regYYstate;

  if (regYYss + regYYstacksize - 1 <= regYYssp)
    {
      /* Get the current used size of the three stacks, in elements.  */
      YYSIZE_T regYYsize = regYYssp - regYYss + 1;

#ifdef regYYoverflow
      {
	/* Give user a chance to reallocate the stack.  Use copies of
	   these so that the &'s don't force the real ones into
	   memory.  */
	YYSTYPE *regYYvs1 = regYYvs;
	regYYtype_int16 *regYYss1 = regYYss;

	/* Each stack pointer address is followed by the size of the
	   data in use in that stack, in bytes.  This used to be a
	   conditional around just the two extra args, but that might
	   be undefined if regYYoverflow is a macro.  */
	regYYoverflow (YY_("memory exhausted"),
		    &regYYss1, regYYsize * sizeof (*regYYssp),
		    &regYYvs1, regYYsize * sizeof (*regYYvsp),
		    &regYYstacksize);

	regYYss = regYYss1;
	regYYvs = regYYvs1;
      }
#else /* no regYYoverflow */
# ifndef YYSTACK_RELOCATE
      goto regYYexhaustedlab;
# else
      /* Extend the stack our own way.  */
      if (YYMAXDEPTH <= regYYstacksize)
	goto regYYexhaustedlab;
      regYYstacksize *= 2;
      if (YYMAXDEPTH < regYYstacksize)
	regYYstacksize = YYMAXDEPTH;

      {
	regYYtype_int16 *regYYss1 = regYYss;
	union regYYalloc *regYYptr =
	  (union regYYalloc *) YYSTACK_ALLOC (YYSTACK_BYTES (regYYstacksize));
	if (! regYYptr)
	  goto regYYexhaustedlab;
	YYSTACK_RELOCATE (regYYss_alloc, regYYss);
	YYSTACK_RELOCATE (regYYvs_alloc, regYYvs);
#  undef YYSTACK_RELOCATE
	if (regYYss1 != regYYssa)
	  YYSTACK_FREE (regYYss1);
      }
# endif
#endif /* no regYYoverflow */

      regYYssp = regYYss + regYYsize - 1;
      regYYvsp = regYYvs + regYYsize - 1;

      YYDPRINTF ((stderr, "Stack size increased to %lu\n",
		  (unsigned long int) regYYstacksize));

      if (regYYss + regYYstacksize - 1 <= regYYssp)
	YYABORT;
    }

  YYDPRINTF ((stderr, "Entering state %d\n", regYYstate));

  if (regYYstate == YYFINAL)
    YYACCEPT;

  goto regYYbackup;

/*-----------.
| regYYbackup.  |
`-----------*/
regYYbackup:

  /* Do appropriate processing given the current state.  Read a
     lookahead token if we need one and don't already have one.  */

  /* First try to decide what to do without reference to lookahead token.  */
  regYYn = regYYpact[regYYstate];
  if (regYYn == YYPACT_NINF)
    goto regYYdefault;

  /* Not known => get a lookahead token if don't already have one.  */

  /* YYCHAR is either YYEMPTY or YYEOF or a valid lookahead symbol.  */
  if (regYYchar == YYEMPTY)
    {
      YYDPRINTF ((stderr, "Reading a token: "));
      regYYchar = YYLEX;
    }

  if (regYYchar <= YYEOF)
    {
      regYYchar = regYYtoken = YYEOF;
      YYDPRINTF ((stderr, "Now at end of input.\n"));
    }
  else
    {
      regYYtoken = YYTRANSLATE (regYYchar);
      YY_SYMBOL_PRINT ("Next token is", regYYtoken, &regYYlval, &regYYlloc);
    }

  /* If the proper action on seeing token YYTOKEN is to reduce or to
     detect an error, take that action.  */
  regYYn += regYYtoken;
  if (regYYn < 0 || YYLAST < regYYn || regYYcheck[regYYn] != regYYtoken)
    goto regYYdefault;
  regYYn = regYYtable[regYYn];
  if (regYYn <= 0)
    {
      if (regYYn == 0 || regYYn == YYTABLE_NINF)
	goto regYYerrlab;
      regYYn = -regYYn;
      goto regYYreduce;
    }

  /* Count tokens shifted since error; after three, turn off error
     status.  */
  if (regYYerrstatus)
    regYYerrstatus--;

  /* Shift the lookahead token.  */
  YY_SYMBOL_PRINT ("Shifting", regYYtoken, &regYYlval, &regYYlloc);

  /* Discard the shifted token.  */
  regYYchar = YYEMPTY;

  regYYstate = regYYn;
  *++regYYvsp = regYYlval;

  goto regYYnewstate;


/*-----------------------------------------------------------.
| regYYdefault -- do the default action for the current state.  |
`-----------------------------------------------------------*/
regYYdefault:
  regYYn = regYYdefact[regYYstate];
  if (regYYn == 0)
    goto regYYerrlab;
  goto regYYreduce;


/*-----------------------------.
| regYYreduce -- Do a reduction.  |
`-----------------------------*/
regYYreduce:
  /* regYYn is the number of a rule to reduce with.  */
  regYYlen = regYYr2[regYYn];

  /* If YYLEN is nonzero, implement the default value of the action:
     `$$ = $1'.

     Otherwise, the following line sets YYVAL to garbage.
     This behavior is undocumented and Bison
     users should not rely upon it.  Assigning to YYVAL
     unconditionally makes the parser a bit smaller, and it avoids a
     GCC warning that YYVAL may be used uninitialized.  */
  regYYval = regYYvsp[1-regYYlen];


  YY_REDUCE_PRINT (regYYn);
  switch (regYYn)
    {
        case 2:

/* Line 1455 of yacc.c  */
#line 45 "regparser.y"
    { (regYYval.my_region) = (regYYvsp[(1) - (3)].my_region); regAddShape( (regYYvsp[(1) - (3)].my_region), regOR, (regYYvsp[(3) - (3)].my_shape) );  }
    break;

  case 3:

/* Line 1455 of yacc.c  */
#line 46 "regparser.y"
    { (regYYval.my_region) = (regYYvsp[(1) - (2)].my_region); regAddShape( (regYYvsp[(1) - (2)].my_region), regOR, (regYYvsp[(2) - (2)].my_shape) );  }
    break;

  case 4:

/* Line 1455 of yacc.c  */
#line 47 "regparser.y"
    { (regYYval.my_region) = (regYYvsp[(1) - (3)].my_region); regAddShape( (regYYvsp[(1) - (3)].my_region), regAND, (regYYvsp[(3) - (3)].my_shape) ); }
    break;

  case 5:

/* Line 1455 of yacc.c  */
#line 48 "regparser.y"
    { (regYYval.my_region) = (regYYvsp[(1) - (3)].my_region); regNegate( (regYYvsp[(3) - (3)].my_shape) ); regAddShape( (regYYvsp[(1) - (3)].my_region), regAND, (regYYvsp[(3) - (3)].my_shape) ); }
    break;

  case 6:

/* Line 1455 of yacc.c  */
#line 49 "regparser.y"
    {
                (regYYval.my_region) = regCreateRegion( NULL, NULL );
                regAddShape( (regYYval.my_region) , regOR, (regYYvsp[(2) - (2)].my_shape) ); 
		my_Gregion = (regYYval.my_region);
             }
    break;

  case 7:

/* Line 1455 of yacc.c  */
#line 54 "regparser.y"
    {
                (regYYval.my_region) = regCreateRegion(NULL, NULL );
                regAddShape( (regYYval.my_region) , regOR, (regYYvsp[(1) - (1)].my_shape) ); 
		my_Gregion = (regYYval.my_region);
             }
    break;

  case 8:

/* Line 1455 of yacc.c  */
#line 59 "regparser.y"
    { regYYerrok; }
    break;

  case 9:

/* Line 1455 of yacc.c  */
#line 66 "regparser.y"
    {  (regYYval.dval) = (regYYvsp[(1) - (1)].dval); world_coord = 0; }
    break;

  case 10:

/* Line 1455 of yacc.c  */
#line 68 "regparser.y"
    {  (regYYval.dval) = (regYYvsp[(1) - (2)].dval); world_coord = 1; }
    break;

  case 11:

/* Line 1455 of yacc.c  */
#line 70 "regparser.y"
    {
           double ndeg, nmin, nsec, nval;
           ndeg = (regYYvsp[(1) - (5)].dval); nmin = (regYYvsp[(3) - (5)].dval); nsec = (regYYvsp[(5) - (5)].dval);
           nval = ( ndeg + nmin / 60.0 + nsec / 3600.0 );
           world_coord = 1;
           (regYYval.dval) = nval * 15.0;
          }
    break;

  case 12:

/* Line 1455 of yacc.c  */
#line 81 "regparser.y"
    {  (regYYval.dval) = (regYYvsp[(1) - (1)].dval); world_coord = 0; }
    break;

  case 13:

/* Line 1455 of yacc.c  */
#line 83 "regparser.y"
    {  (regYYval.dval) = (regYYvsp[(1) - (2)].dval); world_coord = 1; }
    break;

  case 14:

/* Line 1455 of yacc.c  */
#line 85 "regparser.y"
    {
           double ndeg, nmin, nsec, nval;
           ndeg = (regYYvsp[(1) - (5)].dval); nmin = (regYYvsp[(3) - (5)].dval); nsec = (regYYvsp[(5) - (5)].dval);
           nval = ( ndeg + nmin / 60.0 + nsec / 3600.0 );
           world_coord = 1;
           (regYYval.dval) = nval;
          }
    break;

  case 15:

/* Line 1455 of yacc.c  */
#line 93 "regparser.y"
    {  /* Special care in case of -00:00:01 */
           double ndeg, nmin, nsec, nval;
           ndeg = (regYYvsp[(2) - (6)].dval); nmin = (regYYvsp[(4) - (6)].dval); nsec = (regYYvsp[(6) - (6)].dval);
           nval = ( ndeg + nmin / 60.0 + nsec / 3600.0 );
           world_coord = 1;
           (regYYval.dval) = -nval;
          }
    break;

  case 16:

/* Line 1455 of yacc.c  */
#line 104 "regparser.y"
    {  (regYYval.dval) = (regYYvsp[(1) - (1)].dval); }
    break;

  case 17:

/* Line 1455 of yacc.c  */
#line 106 "regparser.y"
    {  (regYYval.dval) = -(regYYvsp[(2) - (2)].dval);}
    break;

  case 18:

/* Line 1455 of yacc.c  */
#line 111 "regparser.y"
    {  (regYYval.dval) = (regYYvsp[(1) - (1)].dval); }
    break;

  case 19:

/* Line 1455 of yacc.c  */
#line 113 "regparser.y"
    {  (regYYval.dval) = (regYYvsp[(1) - (2)].dval) / 60.0; }
    break;

  case 20:

/* Line 1455 of yacc.c  */
#line 115 "regparser.y"
    {  (regYYval.dval) = (regYYvsp[(1) - (2)].dval) / 3600.0; }
    break;

  case 21:

/* Line 1455 of yacc.c  */
#line 117 "regparser.y"
    {  (regYYval.dval) = (regYYvsp[(1) - (2)].dval); }
    break;

  case 22:

/* Line 1455 of yacc.c  */
#line 122 "regparser.y"
    {  (regYYval.dval) = (regYYvsp[(1) - (1)].dval); world_size = 0; }
    break;

  case 23:

/* Line 1455 of yacc.c  */
#line 124 "regparser.y"
    {  (regYYval.dval) = (regYYvsp[(1) - (2)].dval); world_size = 0; }
    break;

  case 24:

/* Line 1455 of yacc.c  */
#line 126 "regparser.y"
    {  (regYYval.dval) = (regYYvsp[(1) - (2)].dval); world_size = 0;  /* This is logical coords; we need to fix this */ }
    break;

  case 25:

/* Line 1455 of yacc.c  */
#line 128 "regparser.y"
    {  (regYYval.dval) = (regYYvsp[(1) - (2)].dval); world_size = 1; }
    break;

  case 26:

/* Line 1455 of yacc.c  */
#line 130 "regparser.y"
    {  (regYYval.dval) = (regYYvsp[(1) - (2)].dval) / 60.0; world_size = 1; }
    break;

  case 27:

/* Line 1455 of yacc.c  */
#line 132 "regparser.y"
    {  (regYYval.dval) = (regYYvsp[(1) - (3)].dval) / 3600.0; world_size = 1; }
    break;

  case 28:

/* Line 1455 of yacc.c  */
#line 134 "regparser.y"
    {  (regYYval.dval) = (regYYvsp[(1) - (2)].dval) / 3600.0; world_size = 1; }
    break;

  case 29:

/* Line 1455 of yacc.c  */
#line 140 "regparser.y"
    { 
       (regYYval.my_shape) = regCreateNewWorldShape( regFIELD, regInclude, NULL, NULL, 0, NULL, NULL, world_coord, world_size );
       if ( (regYYval.my_shape) == NULL ) {
	 my_Gregion = NULL;
	 YYERROR;
       }
     }
    break;

  case 30:

/* Line 1455 of yacc.c  */
#line 148 "regparser.y"
    { 
       (regYYval.my_shape) = regCreateNewWorldShape( regFIELD, regExclude, NULL, NULL, 0, NULL, NULL, world_coord, world_size );
       if ( (regYYval.my_shape) == NULL ) {
	 my_Gregion = NULL;
	 YYERROR;
       }
     }
    break;

  case 31:

/* Line 1455 of yacc.c  */
#line 156 "regparser.y"
    { 
       /* We just ignore the text, treat it as a point */
       double x[1]; double y[1];
       x[0]=(regYYvsp[(3) - (8)].dval); y[0]=(regYYvsp[(5) - (8)].dval);
       (regYYval.my_shape) = regCreateNewWorldShape( regPOINT, regInclude, x, y, 1, NULL, NULL, world_coord, 0 );
       if ( (regYYval.my_shape) == NULL ) {
	 my_Gregion = NULL;
	 YYERROR;
       }
     }
    break;

  case 32:

/* Line 1455 of yacc.c  */
#line 167 "regparser.y"
    { 
       double x[1]; double y[1]; double r[1];
       x[0]=(regYYvsp[(3) - (8)].dval); y[0]=(regYYvsp[(5) - (8)].dval); r[0]=(regYYvsp[(7) - (8)].dval);

       if ( r[0] < 0 ) {
	 fprintf( stderr, "ERROR: circle radius must be positive\n");
	 my_Gregion = NULL;
	 YYERROR;
       } else
	 (regYYval.my_shape) = regCreateNewWorldShape( regCIRCLE, regInclude, x, y, 1, r, NULL, world_coord, world_size );
         if ( (regYYval.my_shape) == NULL ) {
	   my_Gregion = NULL;
	   YYERROR;
         }
     }
    break;

  case 33:

/* Line 1455 of yacc.c  */
#line 183 "regparser.y"
    { 
       double x[1]; double y[1]; double r[1];
       x[0]=(regYYvsp[(4) - (9)].dval); y[0]=(regYYvsp[(6) - (9)].dval); r[0]=(regYYvsp[(8) - (9)].dval);
       if ( r[0] < 0 ) {
	 fprintf( stderr, "ERROR: circle radius must be positive\n");
	 my_Gregion = NULL;
	 YYERROR;
       } else
	 (regYYval.my_shape) = regCreateNewWorldShape( regCIRCLE, regExclude, x, y, 1, r, NULL, world_coord, world_size );
         if ( (regYYval.my_shape) == NULL ) {
  	   my_Gregion = NULL;
	   YYERROR;
         }
     }
    break;

  case 34:

/* Line 1455 of yacc.c  */
#line 198 "regparser.y"
    {
       double x[1]; double y[1]; double r[2];
       x[0]=(regYYvsp[(3) - (10)].dval); y[0]=(regYYvsp[(5) - (10)].dval); r[0]=(regYYvsp[(7) - (10)].dval); r[1]=(regYYvsp[(9) - (10)].dval);

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
	 (regYYval.my_shape) = regCreateNewWorldShape( regANNULUS, regInclude, x, y, 1, r, NULL, world_coord, world_size);
         if ( (regYYval.my_shape) == NULL ) {
 	   my_Gregion = NULL;
	   YYERROR;
         }
     }
    break;

  case 35:

/* Line 1455 of yacc.c  */
#line 219 "regparser.y"
    {
       double x[1]; double y[1]; double r[2];
       x[0]=(regYYvsp[(4) - (11)].dval); y[0]=(regYYvsp[(6) - (11)].dval); r[0]=(regYYvsp[(8) - (11)].dval); r[1]=(regYYvsp[(10) - (11)].dval);

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
	 (regYYval.my_shape) = regCreateNewWorldShape( regANNULUS, regExclude, x, y, 1, r, NULL, world_coord, world_size);
         if ( (regYYval.my_shape) == NULL ) {
 	   my_Gregion = NULL;
	   YYERROR;
         }
     }
    break;

  case 36:

/* Line 1455 of yacc.c  */
#line 240 "regparser.y"
    {
       double x[1]; double y[1]; double r[2];
       x[0]=(regYYvsp[(3) - (10)].dval); y[0]=(regYYvsp[(5) - (10)].dval); r[0]=(regYYvsp[(7) - (10)].dval); r[1]=(regYYvsp[(9) - (10)].dval);

       if (( r[0] < 0 ) || (r[1] < 0 )) {
	 fprintf( stderr, "ERROR: box lengths must be positive\n");
	 my_Gregion = NULL;
	 YYERROR;
       } else
	 (regYYval.my_shape) = regCreateNewWorldShape( regBOX, regInclude, x, y, 1, r, NULL, world_coord, world_size);
         if ( (regYYval.my_shape) == NULL ) {
 	   my_Gregion = NULL;
	   YYERROR;
         }
     }
    break;

  case 37:

/* Line 1455 of yacc.c  */
#line 256 "regparser.y"
    {
       double x[1]; double y[1]; double r[2];
       x[0]=(regYYvsp[(4) - (11)].dval); y[0]=(regYYvsp[(6) - (11)].dval); r[0]=(regYYvsp[(8) - (11)].dval); r[1]=(regYYvsp[(10) - (11)].dval);
       if (( r[0] < 0 ) || (r[1] < 0 )) {
	 fprintf( stderr, "ERROR: box lengths must be positive\n");
	 my_Gregion = NULL;
	 YYERROR;
       } else
	 (regYYval.my_shape) = regCreateNewWorldShape( regBOX, regExclude, x, y, 1, r, NULL, world_coord, world_size);
         if ( (regYYval.my_shape) == NULL ) {
 	   my_Gregion = NULL;
	   YYERROR;
         }
     }
    break;

  case 38:

/* Line 1455 of yacc.c  */
#line 271 "regparser.y"
    {
       double x[1]; double y[1]; double r[2]; double a[1];
       x[0]=(regYYvsp[(3) - (12)].dval); y[0]=(regYYvsp[(5) - (12)].dval); r[0]=(regYYvsp[(7) - (12)].dval); r[1]=(regYYvsp[(9) - (12)].dval); a[0]=(regYYvsp[(11) - (12)].dval);
       if (( r[0] < 0 ) || (r[1] < 0 )) {
	 fprintf( stderr, "ERROR: box lengths must be positive\n");
	 my_Gregion = NULL;
	 YYERROR;
       } else
	 (regYYval.my_shape) = regCreateNewWorldShape( regROTBOX, regInclude, x, y, 1, r, a, world_coord, world_size);
         if ( (regYYval.my_shape) == NULL ) {
 	   my_Gregion = NULL;
	   YYERROR;
         }
     }
    break;

  case 39:

/* Line 1455 of yacc.c  */
#line 286 "regparser.y"
    {
       double x[1]; double y[1]; double r[2]; double a[1];
       x[0]=(regYYvsp[(4) - (13)].dval); y[0]=(regYYvsp[(6) - (13)].dval); r[0]=(regYYvsp[(8) - (13)].dval); r[1]=(regYYvsp[(10) - (13)].dval); a[0]=(regYYvsp[(12) - (13)].dval);
       if (( r[0] < 0 ) || (r[1] < 0 )) {
	 fprintf( stderr, "ERROR: box lengths must be positive\n");
	 my_Gregion = NULL;
	 YYERROR;
       } else
	 (regYYval.my_shape) = regCreateNewWorldShape( regROTBOX, regExclude, x, y, 1, r, a, world_coord, world_size);
         if ( (regYYval.my_shape) == NULL ) {
 	   my_Gregion = NULL;
	   YYERROR;
         }
     }
    break;

  case 40:

/* Line 1455 of yacc.c  */
#line 301 "regparser.y"
    {
       double x[1]; double y[1]; double r[2]; double a[1];
       x[0]=(regYYvsp[(3) - (12)].dval); y[0]=(regYYvsp[(5) - (12)].dval); r[0]=(regYYvsp[(7) - (12)].dval); r[1]=(regYYvsp[(9) - (12)].dval); a[0]=(regYYvsp[(11) - (12)].dval);
       if (( r[0] < 0 ) || (r[1] < 0 )) {
	 fprintf( stderr, "ERROR: box lengths must be positive\n");
	 my_Gregion = NULL;
	 YYERROR;
       } else
	 (regYYval.my_shape) = regCreateNewWorldShape( regROTBOX, regInclude, x, y, 1, r, a, world_coord, world_size);
         if ( (regYYval.my_shape) == NULL ) {
 	   my_Gregion = NULL;
	   YYERROR;
         }
     }
    break;

  case 41:

/* Line 1455 of yacc.c  */
#line 316 "regparser.y"
    {
       double x[1]; double y[1]; double r[2]; double a[1];
       x[0]=(regYYvsp[(4) - (13)].dval); y[0]=(regYYvsp[(6) - (13)].dval); r[0]=(regYYvsp[(8) - (13)].dval); r[1]=(regYYvsp[(10) - (13)].dval); a[0]=(regYYvsp[(12) - (13)].dval);
       if (( r[0] < 0 ) || (r[1] < 0 )) {
	 fprintf( stderr, "ERROR: box lengths must be positive\n");
	 my_Gregion = NULL;
	 YYERROR;
       } else
	 (regYYval.my_shape) = regCreateNewWorldShape( regROTBOX, regExclude, x, y, 1, r, a, world_coord, world_size);
         if ( (regYYval.my_shape) == NULL ) {
 	   my_Gregion = NULL;
	   YYERROR;
         }
     }
    break;

  case 42:

/* Line 1455 of yacc.c  */
#line 331 "regparser.y"
    {
       double x[1]; double y[1]; double a[2]; double r[2];
       x[0]=(regYYvsp[(3) - (14)].dval); y[0]=(regYYvsp[(5) - (14)].dval); r[0]=(regYYvsp[(7) - (14)].dval); r[1]=(regYYvsp[(9) - (14)].dval);  a[0]=(regYYvsp[(11) - (14)].dval); a[1]=(regYYvsp[(13) - (14)].dval);
       if (( r[0] < 0 ) || (r[1] < 0 )) {
	 fprintf( stderr, "ERROR: pie radii must be positive\n");
	 my_Gregion = NULL;
	 YYERROR;
       } else if ( r[0] > r[1] ) {
	 fprintf( stderr, "ERROR: outer annuls radius must be >= inner radius\n");
	 my_Gregion = NULL;
	 YYERROR;
       } else
	 (regYYval.my_shape) = regCreateNewWorldShape( regPIE, regInclude, x, y, 1, r, a, world_coord, world_size );
         if ( (regYYval.my_shape) == NULL ) {
 	   my_Gregion = NULL;
	   YYERROR;
         }
     }
    break;

  case 43:

/* Line 1455 of yacc.c  */
#line 350 "regparser.y"
    {
       double x[1]; double y[1]; double a[2]; double r[2];
       x[0]=(regYYvsp[(4) - (15)].dval); y[0]=(regYYvsp[(6) - (15)].dval); r[0] = (regYYvsp[(8) - (15)].dval); r[1] = (regYYvsp[(10) - (15)].dval); a[0]=(regYYvsp[(12) - (15)].dval); a[1]=(regYYvsp[(14) - (15)].dval);
       if (( r[0] < 0 ) || (r[1] < 0 )) {
	 fprintf( stderr, "ERROR: pie radii must be positive\n");
	 my_Gregion = NULL;
	 YYERROR;
       } else if ( r[0] > r[1] ) {
	 fprintf( stderr, "ERROR: outer annuls radius must be >= inner radius\n");
	 my_Gregion = NULL;
	 YYERROR;
       } else
	 (regYYval.my_shape) = regCreateNewWorldShape( regPIE, regExclude, x, y, 1, r, a, world_coord, world_size );
         if ( (regYYval.my_shape) == NULL ) {
 	   my_Gregion = NULL;
	   YYERROR;
         }
     }
    break;

  case 44:

/* Line 1455 of yacc.c  */
#line 369 "regparser.y"
    {
       double x[1]; double y[1]; double a[2];
       x[0]=(regYYvsp[(3) - (10)].dval); y[0]=(regYYvsp[(5) - (10)].dval); a[0]=(regYYvsp[(7) - (10)].dval); a[1]=(regYYvsp[(9) - (10)].dval);
       (regYYval.my_shape) = regCreateNewWorldShape( regSECTOR, regInclude, x, y, 1, NULL, a, world_coord, 0 );
       if ( (regYYval.my_shape) == NULL ) {
 	 my_Gregion = NULL;
	 YYERROR;
       }
     }
    break;

  case 45:

/* Line 1455 of yacc.c  */
#line 379 "regparser.y"
    {
       double x[1]; double y[1]; double a[2];
       x[0]=(regYYvsp[(4) - (11)].dval); y[0]=(regYYvsp[(6) - (11)].dval); a[0]=(regYYvsp[(8) - (11)].dval); a[1]=(regYYvsp[(10) - (11)].dval);
       (regYYval.my_shape) = regCreateNewWorldShape( regSECTOR, regExclude, x, y, 1, NULL, a, world_coord, 0 );
       if ( (regYYval.my_shape) == NULL ) {
 	 my_Gregion = NULL;
	 YYERROR;
       }
     }
    break;

  case 46:

/* Line 1455 of yacc.c  */
#line 389 "regparser.y"
    {
       double x[1]; double y[1];
       x[0]=(regYYvsp[(3) - (6)].dval); y[0]=(regYYvsp[(5) - (6)].dval);
       (regYYval.my_shape) = regCreateNewWorldShape( regPOINT, regInclude, x, y, 1, NULL, NULL, world_coord, 0 );
       if ( (regYYval.my_shape) == NULL ) {
 	 my_Gregion = NULL;
	 YYERROR;
       }
     }
    break;

  case 47:

/* Line 1455 of yacc.c  */
#line 399 "regparser.y"
    {
       double x[1]; double y[1];
       x[0]=(regYYvsp[(4) - (7)].dval); y[0]=(regYYvsp[(6) - (7)].dval);
       (regYYval.my_shape) = regCreateNewWorldShape( regPOINT, regExclude, x, y, 1, NULL, NULL, world_coord, 0 );
       if ( (regYYval.my_shape) == NULL ) {
 	 my_Gregion = NULL;
	 YYERROR;
       }
     }
    break;

  case 48:

/* Line 1455 of yacc.c  */
#line 409 "regparser.y"
    {
       double x[2]; double y[2];
       x[0]=(regYYvsp[(3) - (10)].dval); x[1]=(regYYvsp[(7) - (10)].dval); y[0]=(regYYvsp[(5) - (10)].dval);y[1]=(regYYvsp[(9) - (10)].dval);
       (regYYval.my_shape) = regCreateNewWorldShape( regRECTANGLE, regInclude, x, y, 2, NULL, NULL, world_coord, 0 );
       if ( (regYYval.my_shape) == NULL ) {
 	 my_Gregion = NULL;
	 YYERROR;
       }
     }
    break;

  case 49:

/* Line 1455 of yacc.c  */
#line 419 "regparser.y"
    {
       double x[2]; double y[2];
       x[0]=(regYYvsp[(4) - (11)].dval); x[1]=(regYYvsp[(8) - (11)].dval); y[0]=(regYYvsp[(6) - (11)].dval); y[1]=(regYYvsp[(10) - (11)].dval);
       (regYYval.my_shape) = regCreateNewWorldShape( regRECTANGLE, regExclude, x, y, 2, NULL, NULL, world_coord, 0 );
       if ( (regYYval.my_shape) == NULL ) {
 	 my_Gregion = NULL;
	 YYERROR;
       }
     }
    break;

  case 50:

/* Line 1455 of yacc.c  */
#line 429 "regparser.y"
    {
       /* RegLine doesn't work correctly; need to calculate angle */
       double x[2]; double y[2];
       x[0]=(regYYvsp[(3) - (10)].dval); x[1]=(regYYvsp[(7) - (10)].dval); y[0]=(regYYvsp[(5) - (10)].dval);y[1]=(regYYvsp[(9) - (10)].dval);
       (regYYval.my_shape) = regCreateNewWorldShape( regRECTANGLE, regInclude, x, y, 2, NULL, NULL, world_coord, 0 );
       if ( (regYYval.my_shape) == NULL ) {
 	 my_Gregion = NULL;
	 YYERROR;
       }
     }
    break;

  case 51:

/* Line 1455 of yacc.c  */
#line 440 "regparser.y"
    {
       double x[2]; double y[2];
       x[0]=(regYYvsp[(4) - (11)].dval); x[1]=(regYYvsp[(8) - (11)].dval); y[0]=(regYYvsp[(6) - (11)].dval); y[1]=(regYYvsp[(10) - (11)].dval);
       (regYYval.my_shape) = regCreateNewWorldShape( regRECTANGLE, regExclude, x, y, 2, NULL, NULL, world_coord, 0 );
       if ( (regYYval.my_shape) == NULL ) {
 	 my_Gregion = NULL;
	 YYERROR;
       }
     }
    break;

  case 52:

/* Line 1455 of yacc.c  */
#line 450 "regparser.y"
    {
       double x[1]; double y[1]; double r[2]; double a[1];
       x[0]=(regYYvsp[(3) - (12)].dval); y[0]=(regYYvsp[(5) - (12)].dval); r[0]=(regYYvsp[(7) - (12)].dval); r[1]=(regYYvsp[(9) - (12)].dval); a[0]=(regYYvsp[(11) - (12)].dval);
       if (( r[0] < 0 ) || (r[1] < 0 )) {
	 fprintf( stderr, "ERROR: ellipse radii must be positive\n");
	 my_Gregion = NULL;
	 YYERROR;
       } else
	 (regYYval.my_shape) = regCreateNewWorldShape( regELLIPSE, regInclude, x, y, 1, r, a, world_coord, world_size );
         if ( (regYYval.my_shape) == NULL ) {
 	   my_Gregion = NULL;
	   YYERROR;
         }
     }
    break;

  case 53:

/* Line 1455 of yacc.c  */
#line 465 "regparser.y"
    {
       double x[1]; double y[1]; double r[2]; double a[1];
       x[0]=(regYYvsp[(4) - (13)].dval); y[0]=(regYYvsp[(6) - (13)].dval); r[0]=(regYYvsp[(8) - (13)].dval); r[1]=(regYYvsp[(10) - (13)].dval); a[0]=(regYYvsp[(12) - (13)].dval);
       if (( r[0] < 0 ) || (r[1] < 0 )) {
	 fprintf( stderr, "ERROR: ellipse radii must be positive\n");
	 my_Gregion = NULL;
	 YYERROR;
       } else
	 (regYYval.my_shape) = regCreateNewWorldShape( regELLIPSE, regExclude, x, y, 1, r, a, world_coord, world_size );
         if ( (regYYval.my_shape) == NULL ) {
 	   my_Gregion = NULL;
	   YYERROR;
         }
     }
    break;

  case 54:

/* Line 1455 of yacc.c  */
#line 480 "regparser.y"
    { 
       (regYYval.my_shape) = regCreateNewWorldShape( regPOLYGON, regInclude, (regYYvsp[(3) - (4)].PolySide).polyX, (regYYvsp[(3) - (4)].PolySide).polyY, 
			    (regYYvsp[(3) - (4)].PolySide).polyS, NULL, NULL, world_coord, 0 );
       if ( (regYYval.my_shape) == NULL ) {
	 my_Gregion = NULL;
	 YYERROR;
       }

       free( (regYYvsp[(3) - (4)].PolySide).polyX );
       free( (regYYvsp[(3) - (4)].PolySide).polyY );
     }
    break;

  case 55:

/* Line 1455 of yacc.c  */
#line 492 "regparser.y"
    { 
       (regYYval.my_shape) = regCreateNewWorldShape( regPOLYGON, regExclude, (regYYvsp[(4) - (5)].PolySide).polyX, (regYYvsp[(4) - (5)].PolySide).polyY, 
			    (regYYvsp[(4) - (5)].PolySide).polyS, NULL, NULL, world_coord, 0 );
       free( (regYYvsp[(4) - (5)].PolySide).polyX );
       free( (regYYvsp[(4) - (5)].PolySide).polyY );
       if ( (regYYval.my_shape) == NULL ) {
	 my_Gregion = NULL;
	 YYERROR;
       }
     }
    break;

  case 56:

/* Line 1455 of yacc.c  */
#line 507 "regparser.y"
    {
     (regYYval.PolySide) = (regYYvsp[(1) - (5)].PolySide);
     (regYYval.PolySide).polyS += 1;
     (regYYval.PolySide).polyX = (double*)realloc( (regYYval.PolySide).polyX, (regYYval.PolySide).polyS *sizeof(double));
     (regYYval.PolySide).polyY = (double*)realloc( (regYYval.PolySide).polyY, (regYYval.PolySide).polyS *sizeof(double));
     (regYYval.PolySide).polyX[(regYYval.PolySide).polyS -1] = (regYYvsp[(3) - (5)].dval);
     (regYYval.PolySide).polyY[(regYYval.PolySide).polyS -1] = (regYYvsp[(5) - (5)].dval);
     }
    break;

  case 57:

/* Line 1455 of yacc.c  */
#line 515 "regparser.y"
    {
     (regYYval.PolySide).polyS = 1;
     (regYYval.PolySide).polyX = (double*)calloc(1,sizeof(double));
     (regYYval.PolySide).polyY = (double*)calloc(1,sizeof(double));
     (regYYval.PolySide).polyX[0] = (regYYvsp[(1) - (3)].dval);
     (regYYval.PolySide).polyY[0] = (regYYvsp[(3) - (3)].dval);
   }
    break;



/* Line 1455 of yacc.c  */
#line 2338 "regparser.tab.c"
      default: break;
    }
  YY_SYMBOL_PRINT ("-> $$ =", regYYr1[regYYn], &regYYval, &regYYloc);

  YYPOPSTACK (regYYlen);
  regYYlen = 0;
  YY_STACK_PRINT (regYYss, regYYssp);

  *++regYYvsp = regYYval;

  /* Now `shift' the result of the reduction.  Determine what state
     that goes to, based on the state we popped back to and the rule
     number reduced by.  */

  regYYn = regYYr1[regYYn];

  regYYstate = regYYpgoto[regYYn - YYNTOKENS] + *regYYssp;
  if (0 <= regYYstate && regYYstate <= YYLAST && regYYcheck[regYYstate] == *regYYssp)
    regYYstate = regYYtable[regYYstate];
  else
    regYYstate = regYYdefgoto[regYYn - YYNTOKENS];

  goto regYYnewstate;


/*------------------------------------.
| regYYerrlab -- here on detecting error |
`------------------------------------*/
regYYerrlab:
  /* If not already recovering from an error, report this error.  */
  if (!regYYerrstatus)
    {
      ++regYYnerrs;
#if ! YYERROR_VERBOSE
      regYYerror (YY_("syntax error"));
#else
      {
	YYSIZE_T regYYsize = regYYsyntax_error (0, regYYstate, regYYchar);
	if (regYYmsg_alloc < regYYsize && regYYmsg_alloc < YYSTACK_ALLOC_MAXIMUM)
	  {
	    YYSIZE_T regYYalloc = 2 * regYYsize;
	    if (! (regYYsize <= regYYalloc && regYYalloc <= YYSTACK_ALLOC_MAXIMUM))
	      regYYalloc = YYSTACK_ALLOC_MAXIMUM;
	    if (regYYmsg != regYYmsgbuf)
	      YYSTACK_FREE (regYYmsg);
	    regYYmsg = (char *) YYSTACK_ALLOC (regYYalloc);
	    if (regYYmsg)
	      regYYmsg_alloc = regYYalloc;
	    else
	      {
		regYYmsg = regYYmsgbuf;
		regYYmsg_alloc = sizeof regYYmsgbuf;
	      }
	  }

	if (0 < regYYsize && regYYsize <= regYYmsg_alloc)
	  {
	    (void) regYYsyntax_error (regYYmsg, regYYstate, regYYchar);
	    regYYerror (regYYmsg);
	  }
	else
	  {
	    regYYerror (YY_("syntax error"));
	    if (regYYsize != 0)
	      goto regYYexhaustedlab;
	  }
      }
#endif
    }



  if (regYYerrstatus == 3)
    {
      /* If just tried and failed to reuse lookahead token after an
	 error, discard it.  */

      if (regYYchar <= YYEOF)
	{
	  /* Return failure if at end of input.  */
	  if (regYYchar == YYEOF)
	    YYABORT;
	}
      else
	{
	  regYYdestruct ("Error: discarding",
		      regYYtoken, &regYYlval);
	  regYYchar = YYEMPTY;
	}
    }

  /* Else will try to reuse lookahead token after shifting the error
     token.  */
  goto regYYerrlab1;


/*---------------------------------------------------.
| regYYerrorlab -- error raised explicitly by YYERROR.  |
`---------------------------------------------------*/
regYYerrorlab:

  /* Pacify compilers like GCC when the user code never invokes
     YYERROR and the label regYYerrorlab therefore never appears in user
     code.  */
  if (/*CONSTCOND*/ 0)
     goto regYYerrorlab;

  /* Do not reclaim the symbols of the rule which action triggered
     this YYERROR.  */
  YYPOPSTACK (regYYlen);
  regYYlen = 0;
  YY_STACK_PRINT (regYYss, regYYssp);
  regYYstate = *regYYssp;
  goto regYYerrlab1;


/*-------------------------------------------------------------.
| regYYerrlab1 -- common code for both syntax error and YYERROR.  |
`-------------------------------------------------------------*/
regYYerrlab1:
  regYYerrstatus = 3;	/* Each real token shifted decrements this.  */

  for (;;)
    {
      regYYn = regYYpact[regYYstate];
      if (regYYn != YYPACT_NINF)
	{
	  regYYn += YYTERROR;
	  if (0 <= regYYn && regYYn <= YYLAST && regYYcheck[regYYn] == YYTERROR)
	    {
	      regYYn = regYYtable[regYYn];
	      if (0 < regYYn)
		break;
	    }
	}

      /* Pop the current state because it cannot handle the error token.  */
      if (regYYssp == regYYss)
	YYABORT;


      regYYdestruct ("Error: popping",
		  regYYstos[regYYstate], regYYvsp);
      YYPOPSTACK (1);
      regYYstate = *regYYssp;
      YY_STACK_PRINT (regYYss, regYYssp);
    }

  *++regYYvsp = regYYlval;


  /* Shift the error token.  */
  YY_SYMBOL_PRINT ("Shifting", regYYstos[regYYn], regYYvsp, regYYlsp);

  regYYstate = regYYn;
  goto regYYnewstate;


/*-------------------------------------.
| regYYacceptlab -- YYACCEPT comes here.  |
`-------------------------------------*/
regYYacceptlab:
  regYYresult = 0;
  goto regYYreturn;

/*-----------------------------------.
| regYYabortlab -- YYABORT comes here.  |
`-----------------------------------*/
regYYabortlab:
  regYYresult = 1;
  goto regYYreturn;

#if !defined(regYYoverflow) || YYERROR_VERBOSE
/*-------------------------------------------------.
| regYYexhaustedlab -- memory exhaustion comes here.  |
`-------------------------------------------------*/
regYYexhaustedlab:
  regYYerror (YY_("memory exhausted"));
  regYYresult = 2;
  /* Fall through.  */
#endif

regYYreturn:
  if (regYYchar != YYEMPTY)
     regYYdestruct ("Cleanup: discarding lookahead",
		 regYYtoken, &regYYlval);
  /* Do not reclaim the symbols of the rule which action triggered
     this YYABORT or YYACCEPT.  */
  YYPOPSTACK (regYYlen);
  YY_STACK_PRINT (regYYss, regYYssp);
  while (regYYssp != regYYss)
    {
      regYYdestruct ("Cleanup: popping",
		  regYYstos[*regYYssp], regYYvsp);
      YYPOPSTACK (1);
    }
#ifndef regYYoverflow
  if (regYYss != regYYssa)
    YYSTACK_FREE (regYYss);
#endif
#if YYERROR_VERBOSE
  if (regYYmsg != regYYmsgbuf)
    YYSTACK_FREE (regYYmsg);
#endif
  /* Make sure YYID is used.  */
  return YYID (regYYresult);
}



/* Line 1675 of yacc.c  */
#line 523 "regparser.y"



regRegion* regParse( char* buf )
{
  regRegion* regptr;
  char* ptr;
  /* my_Gregion is declared externally*/
  my_Gregion = NULL;
  regParseStr = buf;

  if ( !buf )
   return NULL;
  ptr = buf;
  while( *ptr == ' ' || *ptr == '(' ) ptr++;
  if ( !isalpha( *ptr ) && *ptr != '!' )
  {
   /* Not a region ( always begins with alpha or ! ) */
   return NULL;
  }

  regParseStrEnd = buf + strlen( buf );
/*  regYYdebug = 1; */
    regYYparse();
 regptr = my_Gregion;
 return regptr;
}


void regYYerror( char* msg )
{
  
  my_Gregion = NULL;
  return;
}

int regYYwrap( void )
{
  return 1;
}




