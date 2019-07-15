/* A Bison parser, made by GNU Bison 3.0.4.  */

/* Bison implementation for Yacc-like parsers in C

   Copyright (C) 1984, 1989-1990, 2000-2015 Free Software Foundation, Inc.

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
#define YYBISON_VERSION "3.0.4"

/* Skeleton name.  */
#define YYSKELETON_NAME "yacc.c"

/* Pure parsers.  */
#define YYPURE 0

/* Push parsers.  */
#define YYPUSH 0

/* Pull parsers.  */
#define YYPULL 1


/* Substitute the variable and function names.  */
#define yyparse         regYYparse
#define yylex           regYYlex
#define yyerror         regYYerror
#define yydebug         regYYdebug
#define yynerrs         regYYnerrs

#define yylval          regYYlval
#define yychar          regYYchar

/* Copy the first part of user declarations.  */
#line 2 "regparser.y" /* yacc.c:339  */

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

#line 108 "regparser.c" /* yacc.c:339  */

# ifndef YY_NULLPTR
#  if defined __cplusplus && 201103L <= __cplusplus
#   define YY_NULLPTR nullptr
#  else
#   define YY_NULLPTR 0
#  endif
# endif

/* Enabling verbose error messages.  */
#ifdef YYERROR_VERBOSE
# undef YYERROR_VERBOSE
# define YYERROR_VERBOSE 1
#else
# define YYERROR_VERBOSE 0
#endif

/* In a future release of Bison, this section will be replaced
   by #include "y.tab.h".  */
#ifndef YY_REGYY_Y_TAB_H_INCLUDED
# define YY_REGYY_Y_TAB_H_INCLUDED
/* Debug traces.  */
#ifndef YYDEBUG
# define YYDEBUG 0
#endif
#if YYDEBUG
extern int regYYdebug;
#endif

/* Token type.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
  enum yytokentype
  {
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

/* Value type.  */
#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED

union YYSTYPE
{
#line 45 "regparser.y" /* yacc.c:355  */

  double dval;
  char str[1024];
  regRegion* my_region;
  regShape* my_shape;
  struct polyside
  { double *polyX; double *polyY; long polyS; } PolySide;

#line 209 "regparser.c" /* yacc.c:355  */
};

typedef union YYSTYPE YYSTYPE;
# define YYSTYPE_IS_TRIVIAL 1
# define YYSTYPE_IS_DECLARED 1
#endif


extern YYSTYPE regYYlval;

int regYYparse (void);

#endif /* !YY_REGYY_Y_TAB_H_INCLUDED  */

/* Copy the second part of user declarations.  */

#line 226 "regparser.c" /* yacc.c:358  */

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
#else
typedef signed char yytype_int8;
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
# elif ! defined YYSIZE_T
#  include <stddef.h> /* INFRINGES ON USER NAME SPACE */
#  define YYSIZE_T size_t
# else
#  define YYSIZE_T unsigned int
# endif
#endif

#define YYSIZE_MAXIMUM ((YYSIZE_T) -1)

#ifndef YY_
# if defined YYENABLE_NLS && YYENABLE_NLS
#  if ENABLE_NLS
#   include <libintl.h> /* INFRINGES ON USER NAME SPACE */
#   define YY_(Msgid) dgettext ("bison-runtime", Msgid)
#  endif
# endif
# ifndef YY_
#  define YY_(Msgid) Msgid
# endif
#endif

#ifndef YY_ATTRIBUTE
# if (defined __GNUC__                                               \
      && (2 < __GNUC__ || (__GNUC__ == 2 && 96 <= __GNUC_MINOR__)))  \
     || defined __SUNPRO_C && 0x5110 <= __SUNPRO_C
#  define YY_ATTRIBUTE(Spec) __attribute__(Spec)
# else
#  define YY_ATTRIBUTE(Spec) /* empty */
# endif
#endif

#ifndef YY_ATTRIBUTE_PURE
# define YY_ATTRIBUTE_PURE   YY_ATTRIBUTE ((__pure__))
#endif

#ifndef YY_ATTRIBUTE_UNUSED
# define YY_ATTRIBUTE_UNUSED YY_ATTRIBUTE ((__unused__))
#endif

#if !defined _Noreturn \
     && (!defined __STDC_VERSION__ || __STDC_VERSION__ < 201112)
# if defined _MSC_VER && 1200 <= _MSC_VER
#  define _Noreturn __declspec (noreturn)
# else
#  define _Noreturn YY_ATTRIBUTE ((__noreturn__))
# endif
#endif

/* Suppress unused-variable warnings by "using" E.  */
#if ! defined lint || defined __GNUC__
# define YYUSE(E) ((void) (E))
#else
# define YYUSE(E) /* empty */
#endif

#if defined __GNUC__ && 407 <= __GNUC__ * 100 + __GNUC_MINOR__
/* Suppress an incorrect diagnostic about yylval being uninitialized.  */
# define YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN \
    _Pragma ("GCC diagnostic push") \
    _Pragma ("GCC diagnostic ignored \"-Wuninitialized\"")\
    _Pragma ("GCC diagnostic ignored \"-Wmaybe-uninitialized\"")
# define YY_IGNORE_MAYBE_UNINITIALIZED_END \
    _Pragma ("GCC diagnostic pop")
#else
# define YY_INITIAL_VALUE(Value) Value
#endif
#ifndef YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
# define YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
# define YY_IGNORE_MAYBE_UNINITIALIZED_END
#endif
#ifndef YY_INITIAL_VALUE
# define YY_INITIAL_VALUE(Value) /* Nothing. */
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
#    if ! defined _ALLOCA_H && ! defined EXIT_SUCCESS
#     include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
      /* Use EXIT_SUCCESS as a witness for stdlib.h.  */
#     ifndef EXIT_SUCCESS
#      define EXIT_SUCCESS 0
#     endif
#    endif
#   endif
#  endif
# endif

# ifdef YYSTACK_ALLOC
   /* Pacify GCC's 'empty if-body' warning.  */
#  define YYSTACK_FREE(Ptr) do { /* empty */; } while (0)
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
#  if (defined __cplusplus && ! defined EXIT_SUCCESS \
       && ! ((defined YYMALLOC || defined malloc) \
             && (defined YYFREE || defined free)))
#   include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#   ifndef EXIT_SUCCESS
#    define EXIT_SUCCESS 0
#   endif
#  endif
#  ifndef YYMALLOC
#   define YYMALLOC malloc
#   if ! defined malloc && ! defined EXIT_SUCCESS
void *malloc (YYSIZE_T); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
#  ifndef YYFREE
#   define YYFREE free
#   if ! defined free && ! defined EXIT_SUCCESS
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

# define YYCOPY_NEEDED 1

/* Relocate STACK from its old location to the new one.  The
   local variables YYSIZE and YYSTACKSIZE give the old and new number of
   elements in the stack, and YYPTR gives the new location of the
   stack.  Advance YYPTR to a properly aligned location for the next
   stack.  */
# define YYSTACK_RELOCATE(Stack_alloc, Stack)                           \
    do                                                                  \
      {                                                                 \
        YYSIZE_T yynewbytes;                                            \
        YYCOPY (&yyptr->Stack_alloc, Stack, yysize);                    \
        Stack = &yyptr->Stack_alloc;                                    \
        yynewbytes = yystacksize * sizeof (*Stack) + YYSTACK_GAP_MAXIMUM; \
        yyptr += yynewbytes / sizeof (*yyptr);                          \
      }                                                                 \
    while (0)

#endif

#if defined YYCOPY_NEEDED && YYCOPY_NEEDED
/* Copy COUNT objects from SRC to DST.  The source and destination do
   not overlap.  */
# ifndef YYCOPY
#  if defined __GNUC__ && 1 < __GNUC__
#   define YYCOPY(Dst, Src, Count) \
      __builtin_memcpy (Dst, Src, (Count) * sizeof (*(Src)))
#  else
#   define YYCOPY(Dst, Src, Count)              \
      do                                        \
        {                                       \
          YYSIZE_T yyi;                         \
          for (yyi = 0; yyi < (Count); yyi++)   \
            (Dst)[yyi] = (Src)[yyi];            \
        }                                       \
      while (0)
#  endif
# endif
#endif /* !YYCOPY_NEEDED */

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
/* YYNSTATES -- Number of states.  */
#define YYNSTATES  358

/* YYTRANSLATE[YYX] -- Symbol number corresponding to YYX as returned
   by yylex, with out-of-bounds checking.  */
#define YYUNDEFTOK  2
#define YYMAXUTOK   282

#define YYTRANSLATE(YYX)                                                \
  ((unsigned int) (YYX) <= YYMAXUTOK ? yytranslate[YYX] : YYUNDEFTOK)

/* YYTRANSLATE[TOKEN-NUM] -- Symbol number corresponding to TOKEN-NUM
   as returned by yylex, without out-of-bounds checking.  */
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
  /* YYRLINE[YYN] -- Source line where rule number YYN was defined.  */
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

#if YYDEBUG || YYERROR_VERBOSE || 0
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
  "reg_value", "reg_angle", "reg_size", "reg_shape", "reg_ordered_pair", YY_NULLPTR
};
#endif

# ifdef YYPRINT
/* YYTOKNUM[NUM] -- (External) token number corresponding to the
   (internal) symbol number NUM (which must be that of a token).  */
static const yytype_uint16 yytoknum[] =
{
       0,   256,   257,   258,   259,   260,   261,   262,   263,   264,
     265,   266,   267,   268,   269,   270,   271,   272,   273,   274,
     275,   276,   277,   278,   279,   280,   281,   282,   100,    58,
      39,    34,   112,   105,    40,    41,    44
};
# endif

#define YYPACT_NINF -32

#define yypact_value_is_default(Yystate) \
  (!!((Yystate) == (-32)))

#define YYTABLE_NINF -1

#define yytable_value_is_error(Yytable_value) \
  0

  /* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
     STATE-NUM.  */
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

  /* YYDEFACT[STATE-NUM] -- Default reduction number in state STATE-NUM.
     Performed when YYTABLE does not specify something else to do.  Zero
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

  /* YYPGOTO[NTERM-NUM].  */
static const yytype_int16 yypgoto[] =
{
     -32,   -32,   150,   228,   -23,   -31,   197,    11,    62
};

  /* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int16 yydefgoto[] =
{
      -1,    20,    64,   133,   198,   208,   199,    21,    65
};

  /* YYTABLE[YYPACT[STATE-NUM]] -- What to do in state STATE-NUM.  If
     positive, shift that token.  If negative, reduce the rule whose
     number is the opposite.  If YYTABLE_NINF, syntax error.  */
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

  /* YYR2[YYN] -- Number of symbols on the right hand side of rule YYN.  */
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


#define yyerrok         (yyerrstatus = 0)
#define yyclearin       (yychar = YYEMPTY)
#define YYEMPTY         (-2)
#define YYEOF           0

#define YYACCEPT        goto yyacceptlab
#define YYABORT         goto yyabortlab
#define YYERROR         goto yyerrorlab


#define YYRECOVERING()  (!!yyerrstatus)

#define YYBACKUP(Token, Value)                                  \
do                                                              \
  if (yychar == YYEMPTY)                                        \
    {                                                           \
      yychar = (Token);                                         \
      yylval = (Value);                                         \
      YYPOPSTACK (yylen);                                       \
      yystate = *yyssp;                                         \
      goto yybackup;                                            \
    }                                                           \
  else                                                          \
    {                                                           \
      yyerror (YY_("syntax error: cannot back up")); \
      YYERROR;                                                  \
    }                                                           \
while (0)

/* Error token number */
#define YYTERROR        1
#define YYERRCODE       256



/* Enable debugging if requested.  */
#if YYDEBUG

# ifndef YYFPRINTF
#  include <stdio.h> /* INFRINGES ON USER NAME SPACE */
#  define YYFPRINTF fprintf
# endif

# define YYDPRINTF(Args)                        \
do {                                            \
  if (yydebug)                                  \
    YYFPRINTF Args;                             \
} while (0)

/* This macro is provided for backward compatibility. */
#ifndef YY_LOCATION_PRINT
# define YY_LOCATION_PRINT(File, Loc) ((void) 0)
#endif


# define YY_SYMBOL_PRINT(Title, Type, Value, Location)                    \
do {                                                                      \
  if (yydebug)                                                            \
    {                                                                     \
      YYFPRINTF (stderr, "%s ", Title);                                   \
      yy_symbol_print (stderr,                                            \
                  Type, Value); \
      YYFPRINTF (stderr, "\n");                                           \
    }                                                                     \
} while (0)


/*----------------------------------------.
| Print this symbol's value on YYOUTPUT.  |
`----------------------------------------*/

static void
yy_symbol_value_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep)
{
  FILE *yyo = yyoutput;
  YYUSE (yyo);
  if (!yyvaluep)
    return;
# ifdef YYPRINT
  if (yytype < YYNTOKENS)
    YYPRINT (yyoutput, yytoknum[yytype], *yyvaluep);
# endif
  YYUSE (yytype);
}


/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
`--------------------------------*/

static void
yy_symbol_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep)
{
  YYFPRINTF (yyoutput, "%s %s (",
             yytype < YYNTOKENS ? "token" : "nterm", yytname[yytype]);

  yy_symbol_value_print (yyoutput, yytype, yyvaluep);
  YYFPRINTF (yyoutput, ")");
}

/*------------------------------------------------------------------.
| yy_stack_print -- Print the state stack from its BOTTOM up to its |
| TOP (included).                                                   |
`------------------------------------------------------------------*/

static void
yy_stack_print (yytype_int16 *yybottom, yytype_int16 *yytop)
{
  YYFPRINTF (stderr, "Stack now");
  for (; yybottom <= yytop; yybottom++)
    {
      int yybot = *yybottom;
      YYFPRINTF (stderr, " %d", yybot);
    }
  YYFPRINTF (stderr, "\n");
}

# define YY_STACK_PRINT(Bottom, Top)                            \
do {                                                            \
  if (yydebug)                                                  \
    yy_stack_print ((Bottom), (Top));                           \
} while (0)


/*------------------------------------------------.
| Report that the YYRULE is going to be reduced.  |
`------------------------------------------------*/

static void
yy_reduce_print (yytype_int16 *yyssp, YYSTYPE *yyvsp, int yyrule)
{
  unsigned long int yylno = yyrline[yyrule];
  int yynrhs = yyr2[yyrule];
  int yyi;
  YYFPRINTF (stderr, "Reducing stack by rule %d (line %lu):\n",
             yyrule - 1, yylno);
  /* The symbols being reduced.  */
  for (yyi = 0; yyi < yynrhs; yyi++)
    {
      YYFPRINTF (stderr, "   $%d = ", yyi + 1);
      yy_symbol_print (stderr,
                       yystos[yyssp[yyi + 1 - yynrhs]],
                       &(yyvsp[(yyi + 1) - (yynrhs)])
                                              );
      YYFPRINTF (stderr, "\n");
    }
}

# define YY_REDUCE_PRINT(Rule)          \
do {                                    \
  if (yydebug)                          \
    yy_reduce_print (yyssp, yyvsp, Rule); \
} while (0)

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
#ifndef YYINITDEPTH
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
static YYSIZE_T
yystrlen (const char *yystr)
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
static char *
yystpcpy (char *yydest, const char *yysrc)
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

/* Copy into *YYMSG, which is of size *YYMSG_ALLOC, an error message
   about the unexpected token YYTOKEN for the state stack whose top is
   YYSSP.

   Return 0 if *YYMSG was successfully written.  Return 1 if *YYMSG is
   not large enough to hold the message.  In that case, also set
   *YYMSG_ALLOC to the required number of bytes.  Return 2 if the
   required number of bytes is too large to store.  */
static int
yysyntax_error (YYSIZE_T *yymsg_alloc, char **yymsg,
                yytype_int16 *yyssp, int yytoken)
{
  YYSIZE_T yysize0 = yytnamerr (YY_NULLPTR, yytname[yytoken]);
  YYSIZE_T yysize = yysize0;
  enum { YYERROR_VERBOSE_ARGS_MAXIMUM = 5 };
  /* Internationalized format string. */
  const char *yyformat = YY_NULLPTR;
  /* Arguments of yyformat. */
  char const *yyarg[YYERROR_VERBOSE_ARGS_MAXIMUM];
  /* Number of reported tokens (one for the "unexpected", one per
     "expected"). */
  int yycount = 0;

  /* There are many possibilities here to consider:
     - If this state is a consistent state with a default action, then
       the only way this function was invoked is if the default action
       is an error action.  In that case, don't check for expected
       tokens because there are none.
     - The only way there can be no lookahead present (in yychar) is if
       this state is a consistent state with a default action.  Thus,
       detecting the absence of a lookahead is sufficient to determine
       that there is no unexpected or expected token to report.  In that
       case, just report a simple "syntax error".
     - Don't assume there isn't a lookahead just because this state is a
       consistent state with a default action.  There might have been a
       previous inconsistent state, consistent state with a non-default
       action, or user semantic action that manipulated yychar.
     - Of course, the expected token list depends on states to have
       correct lookahead information, and it depends on the parser not
       to perform extra reductions after fetching a lookahead from the
       scanner and before detecting a syntax error.  Thus, state merging
       (from LALR or IELR) and default reductions corrupt the expected
       token list.  However, the list is correct for canonical LR with
       one exception: it will still contain any token that will not be
       accepted due to an error action in a later state.
  */
  if (yytoken != YYEMPTY)
    {
      int yyn = yypact[*yyssp];
      yyarg[yycount++] = yytname[yytoken];
      if (!yypact_value_is_default (yyn))
        {
          /* Start YYX at -YYN if negative to avoid negative indexes in
             YYCHECK.  In other words, skip the first -YYN actions for
             this state because they are default actions.  */
          int yyxbegin = yyn < 0 ? -yyn : 0;
          /* Stay within bounds of both yycheck and yytname.  */
          int yychecklim = YYLAST - yyn + 1;
          int yyxend = yychecklim < YYNTOKENS ? yychecklim : YYNTOKENS;
          int yyx;

          for (yyx = yyxbegin; yyx < yyxend; ++yyx)
            if (yycheck[yyx + yyn] == yyx && yyx != YYTERROR
                && !yytable_value_is_error (yytable[yyx + yyn]))
              {
                if (yycount == YYERROR_VERBOSE_ARGS_MAXIMUM)
                  {
                    yycount = 1;
                    yysize = yysize0;
                    break;
                  }
                yyarg[yycount++] = yytname[yyx];
                {
                  YYSIZE_T yysize1 = yysize + yytnamerr (YY_NULLPTR, yytname[yyx]);
                  if (! (yysize <= yysize1
                         && yysize1 <= YYSTACK_ALLOC_MAXIMUM))
                    return 2;
                  yysize = yysize1;
                }
              }
        }
    }

  switch (yycount)
    {
# define YYCASE_(N, S)                      \
      case N:                               \
        yyformat = S;                       \
      break
      YYCASE_(0, YY_("syntax error"));
      YYCASE_(1, YY_("syntax error, unexpected %s"));
      YYCASE_(2, YY_("syntax error, unexpected %s, expecting %s"));
      YYCASE_(3, YY_("syntax error, unexpected %s, expecting %s or %s"));
      YYCASE_(4, YY_("syntax error, unexpected %s, expecting %s or %s or %s"));
      YYCASE_(5, YY_("syntax error, unexpected %s, expecting %s or %s or %s or %s"));
# undef YYCASE_
    }

  {
    YYSIZE_T yysize1 = yysize + yystrlen (yyformat);
    if (! (yysize <= yysize1 && yysize1 <= YYSTACK_ALLOC_MAXIMUM))
      return 2;
    yysize = yysize1;
  }

  if (*yymsg_alloc < yysize)
    {
      *yymsg_alloc = 2 * yysize;
      if (! (yysize <= *yymsg_alloc
             && *yymsg_alloc <= YYSTACK_ALLOC_MAXIMUM))
        *yymsg_alloc = YYSTACK_ALLOC_MAXIMUM;
      return 1;
    }

  /* Avoid sprintf, as that infringes on the user's name space.
     Don't have undefined behavior even if the translation
     produced a string with the wrong number of "%s"s.  */
  {
    char *yyp = *yymsg;
    int yyi = 0;
    while ((*yyp = *yyformat) != '\0')
      if (*yyp == '%' && yyformat[1] == 's' && yyi < yycount)
        {
          yyp += yytnamerr (yyp, yyarg[yyi++]);
          yyformat += 2;
        }
      else
        {
          yyp++;
          yyformat++;
        }
  }
  return 0;
}
#endif /* YYERROR_VERBOSE */

/*-----------------------------------------------.
| Release the memory associated to this symbol.  |
`-----------------------------------------------*/

static void
yydestruct (const char *yymsg, int yytype, YYSTYPE *yyvaluep)
{
  YYUSE (yyvaluep);
  if (!yymsg)
    yymsg = "Deleting";
  YY_SYMBOL_PRINT (yymsg, yytype, yyvaluep, yylocationp);

  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  YYUSE (yytype);
  YY_IGNORE_MAYBE_UNINITIALIZED_END
}




/* The lookahead symbol.  */
int yychar;

/* The semantic value of the lookahead symbol.  */
YYSTYPE yylval;
/* Number of syntax errors so far.  */
int yynerrs;


/*----------.
| yyparse.  |
`----------*/

int
yyparse (void)
{
    int yystate;
    /* Number of tokens to shift before error messages enabled.  */
    int yyerrstatus;

    /* The stacks and their tools:
       'yyss': related to states.
       'yyvs': related to semantic values.

       Refer to the stacks through separate pointers, to allow yyoverflow
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
  int yytoken = 0;
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

  yyssp = yyss = yyssa;
  yyvsp = yyvs = yyvsa;
  yystacksize = YYINITDEPTH;

  YYDPRINTF ((stderr, "Starting parse\n"));

  yystate = 0;
  yyerrstatus = 0;
  yynerrs = 0;
  yychar = YYEMPTY; /* Cause a token to be read.  */
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
  if (yypact_value_is_default (yyn))
    goto yydefault;

  /* Not known => get a lookahead token if don't already have one.  */

  /* YYCHAR is either YYEMPTY or YYEOF or a valid lookahead symbol.  */
  if (yychar == YYEMPTY)
    {
      YYDPRINTF ((stderr, "Reading a token: "));
      yychar = yylex ();
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
      if (yytable_value_is_error (yyn))
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
  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  *++yyvsp = yylval;
  YY_IGNORE_MAYBE_UNINITIALIZED_END

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
     '$$ = $1'.

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
#line 66 "regparser.y" /* yacc.c:1646  */
    { (yyval.my_region) = (yyvsp[-2].my_region); regAddShape( (yyvsp[-2].my_region), regOR, (yyvsp[0].my_shape) );  }
#line 1542 "regparser.c" /* yacc.c:1646  */
    break;

  case 3:
#line 67 "regparser.y" /* yacc.c:1646  */
    { (yyval.my_region) = (yyvsp[-1].my_region); regAddShape( (yyvsp[-1].my_region), regOR, (yyvsp[0].my_shape) );  }
#line 1548 "regparser.c" /* yacc.c:1646  */
    break;

  case 4:
#line 68 "regparser.y" /* yacc.c:1646  */
    { (yyval.my_region) = (yyvsp[-2].my_region); regAddShape( (yyvsp[-2].my_region), regAND, (yyvsp[0].my_shape) ); }
#line 1554 "regparser.c" /* yacc.c:1646  */
    break;

  case 5:
#line 69 "regparser.y" /* yacc.c:1646  */
    { (yyval.my_region) = (yyvsp[-2].my_region); regNegate( (yyvsp[0].my_shape) ); regAddShape( (yyvsp[-2].my_region), regAND, (yyvsp[0].my_shape) ); }
#line 1560 "regparser.c" /* yacc.c:1646  */
    break;

  case 6:
#line 70 "regparser.y" /* yacc.c:1646  */
    {
                (yyval.my_region) = regCreateRegion( NULL, NULL );
                regAddShape( (yyval.my_region) , regOR, (yyvsp[0].my_shape) ); 
		my_Gregion = (yyval.my_region);
             }
#line 1570 "regparser.c" /* yacc.c:1646  */
    break;

  case 7:
#line 75 "regparser.y" /* yacc.c:1646  */
    {
                (yyval.my_region) = regCreateRegion(NULL, NULL );
                regAddShape( (yyval.my_region) , regOR, (yyvsp[0].my_shape) ); 
		        my_Gregion = (yyval.my_region);
             }
#line 1580 "regparser.c" /* yacc.c:1646  */
    break;

  case 8:
#line 80 "regparser.y" /* yacc.c:1646  */
    { yyerrok; }
#line 1586 "regparser.c" /* yacc.c:1646  */
    break;

  case 9:
#line 87 "regparser.y" /* yacc.c:1646  */
    {  (yyval.dval) = (yyvsp[0].dval); world_coord = RC_UNK; }
#line 1592 "regparser.c" /* yacc.c:1646  */
    break;

  case 10:
#line 89 "regparser.y" /* yacc.c:1646  */
    {  (yyval.dval) = (yyvsp[-1].dval); world_coord = RC_WORLD; }
#line 1598 "regparser.c" /* yacc.c:1646  */
    break;

  case 11:
#line 91 "regparser.y" /* yacc.c:1646  */
    {
           double ndeg, nmin, nsec, nval;
           ndeg = (yyvsp[-4].dval); nmin = (yyvsp[-2].dval); nsec = (yyvsp[0].dval);
           nval = ( ndeg + nmin / 60.0 + nsec / 3600.0 );
           world_coord = RC_WORLD;
           (yyval.dval) = nval * 15.0;
          }
#line 1610 "regparser.c" /* yacc.c:1646  */
    break;

  case 12:
#line 102 "regparser.y" /* yacc.c:1646  */
    {  (yyval.dval) = (yyvsp[0].dval); world_coord = RC_UNK; }
#line 1616 "regparser.c" /* yacc.c:1646  */
    break;

  case 13:
#line 104 "regparser.y" /* yacc.c:1646  */
    {  (yyval.dval) = (yyvsp[-1].dval); world_coord = RC_WORLD; }
#line 1622 "regparser.c" /* yacc.c:1646  */
    break;

  case 14:
#line 106 "regparser.y" /* yacc.c:1646  */
    {
           double ndeg, nmin, nsec, nval;
           ndeg = (yyvsp[-4].dval); nmin = (yyvsp[-2].dval); nsec = (yyvsp[0].dval);
           nval = ( ndeg + nmin / 60.0 + nsec / 3600.0 );
           world_coord = RC_WORLD;
           (yyval.dval) = nval;
          }
#line 1634 "regparser.c" /* yacc.c:1646  */
    break;

  case 15:
#line 114 "regparser.y" /* yacc.c:1646  */
    {  /* Special care in case of -00:00:01 */
           double ndeg, nmin, nsec, nval;
           ndeg = (yyvsp[-4].dval); nmin = (yyvsp[-2].dval); nsec = (yyvsp[0].dval);
           nval = ( ndeg + nmin / 60.0 + nsec / 3600.0 );
           world_coord = RC_WORLD;
           (yyval.dval) = -nval;
          }
#line 1646 "regparser.c" /* yacc.c:1646  */
    break;

  case 16:
#line 125 "regparser.y" /* yacc.c:1646  */
    {  (yyval.dval) = (yyvsp[0].dval); }
#line 1652 "regparser.c" /* yacc.c:1646  */
    break;

  case 17:
#line 127 "regparser.y" /* yacc.c:1646  */
    {  (yyval.dval) = -(yyvsp[0].dval);}
#line 1658 "regparser.c" /* yacc.c:1646  */
    break;

  case 18:
#line 132 "regparser.y" /* yacc.c:1646  */
    {  (yyval.dval) = (yyvsp[0].dval); }
#line 1664 "regparser.c" /* yacc.c:1646  */
    break;

  case 19:
#line 134 "regparser.y" /* yacc.c:1646  */
    {  (yyval.dval) = (yyvsp[-1].dval) / 60.0; }
#line 1670 "regparser.c" /* yacc.c:1646  */
    break;

  case 20:
#line 136 "regparser.y" /* yacc.c:1646  */
    {  (yyval.dval) = (yyvsp[-1].dval) / 3600.0; }
#line 1676 "regparser.c" /* yacc.c:1646  */
    break;

  case 21:
#line 138 "regparser.y" /* yacc.c:1646  */
    {  (yyval.dval) = (yyvsp[-1].dval); }
#line 1682 "regparser.c" /* yacc.c:1646  */
    break;

  case 22:
#line 143 "regparser.y" /* yacc.c:1646  */
    {  (yyval.dval) = (yyvsp[0].dval); world_size = RC_UNK; }
#line 1688 "regparser.c" /* yacc.c:1646  */
    break;

  case 23:
#line 145 "regparser.y" /* yacc.c:1646  */
    {  (yyval.dval) = (yyvsp[-1].dval); world_size = RC_PHYSICAL; }
#line 1694 "regparser.c" /* yacc.c:1646  */
    break;

  case 24:
#line 147 "regparser.y" /* yacc.c:1646  */
    {  (yyval.dval) = (yyvsp[-1].dval); world_size = RC_LOGICAL; }
#line 1700 "regparser.c" /* yacc.c:1646  */
    break;

  case 25:
#line 149 "regparser.y" /* yacc.c:1646  */
    {  (yyval.dval) = (yyvsp[-1].dval); world_size = RC_WORLD; }
#line 1706 "regparser.c" /* yacc.c:1646  */
    break;

  case 26:
#line 151 "regparser.y" /* yacc.c:1646  */
    {  (yyval.dval) = (yyvsp[-1].dval) / 60.0; world_size = RC_WORLD; }
#line 1712 "regparser.c" /* yacc.c:1646  */
    break;

  case 27:
#line 153 "regparser.y" /* yacc.c:1646  */
    {  (yyval.dval) = (yyvsp[-2].dval) / 3600.0; world_size = RC_WORLD; }
#line 1718 "regparser.c" /* yacc.c:1646  */
    break;

  case 28:
#line 155 "regparser.y" /* yacc.c:1646  */
    {  (yyval.dval) = (yyvsp[-1].dval) / 3600.0; world_size = RC_WORLD; }
#line 1724 "regparser.c" /* yacc.c:1646  */
    break;

  case 29:
#line 162 "regparser.y" /* yacc.c:1646  */
    { 
       (yyval.my_shape) = regCreateNewWorldShape( regFIELD, regInclude, NULL, NULL, 0, NULL, NULL, world_coord, world_size );
       if ( (yyval.my_shape) == NULL ) {
	 my_Gregion = NULL;
	 YYERROR;
       }
     }
#line 1736 "regparser.c" /* yacc.c:1646  */
    break;

  case 30:
#line 170 "regparser.y" /* yacc.c:1646  */
    { 
       (yyval.my_shape) = regCreateNewWorldShape( regFIELD, regExclude, NULL, NULL, 0, NULL, NULL, world_coord, world_size );
       if ( (yyval.my_shape) == NULL ) {
	 my_Gregion = NULL;
	 YYERROR;
       }
     }
#line 1748 "regparser.c" /* yacc.c:1646  */
    break;

  case 31:
#line 178 "regparser.y" /* yacc.c:1646  */
    { 
       /* We just ignore the text, treat it as a point */
       double x[1]; double y[1];
       x[0]=(yyvsp[-5].dval); y[0]=(yyvsp[-3].dval);
       (yyval.my_shape) = regCreateNewWorldShape( regPOINT, regInclude, x, y, 1, NULL, NULL, world_coord, 0 );
       if ( (yyval.my_shape) == NULL ) {
	 my_Gregion = NULL;
	 YYERROR;
       }
     }
#line 1763 "regparser.c" /* yacc.c:1646  */
    break;

  case 32:
#line 189 "regparser.y" /* yacc.c:1646  */
    { 
       double x[1]; double y[1]; double r[1];
       x[0]=(yyvsp[-5].dval); y[0]=(yyvsp[-3].dval); r[0]=(yyvsp[-1].dval);

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
#line 1783 "regparser.c" /* yacc.c:1646  */
    break;

  case 33:
#line 205 "regparser.y" /* yacc.c:1646  */
    { 
       double x[1]; double y[1]; double r[1];
       x[0]=(yyvsp[-5].dval); y[0]=(yyvsp[-3].dval); r[0]=(yyvsp[-1].dval);
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
#line 1802 "regparser.c" /* yacc.c:1646  */
    break;

  case 34:
#line 220 "regparser.y" /* yacc.c:1646  */
    {
       double x[1]; double y[1]; double r[2];
       x[0]=(yyvsp[-7].dval); y[0]=(yyvsp[-5].dval); r[0]=(yyvsp[-3].dval); r[1]=(yyvsp[-1].dval);

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
#line 1827 "regparser.c" /* yacc.c:1646  */
    break;

  case 35:
#line 241 "regparser.y" /* yacc.c:1646  */
    {
       double x[1]; double y[1]; double r[2];
       x[0]=(yyvsp[-7].dval); y[0]=(yyvsp[-5].dval); r[0]=(yyvsp[-3].dval); r[1]=(yyvsp[-1].dval);

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
#line 1852 "regparser.c" /* yacc.c:1646  */
    break;

  case 36:
#line 262 "regparser.y" /* yacc.c:1646  */
    {
       double x[1]; double y[1]; double r[2];
       x[0]=(yyvsp[-7].dval); y[0]=(yyvsp[-5].dval); r[0]=(yyvsp[-3].dval); r[1]=(yyvsp[-1].dval);

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
#line 1872 "regparser.c" /* yacc.c:1646  */
    break;

  case 37:
#line 278 "regparser.y" /* yacc.c:1646  */
    {
       double x[1]; double y[1]; double r[2];
       x[0]=(yyvsp[-7].dval); y[0]=(yyvsp[-5].dval); r[0]=(yyvsp[-3].dval); r[1]=(yyvsp[-1].dval);
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
#line 1891 "regparser.c" /* yacc.c:1646  */
    break;

  case 38:
#line 293 "regparser.y" /* yacc.c:1646  */
    {
       double x[1]; double y[1]; double r[2]; double a[1];
       x[0]=(yyvsp[-9].dval); y[0]=(yyvsp[-7].dval); r[0]=(yyvsp[-5].dval); r[1]=(yyvsp[-3].dval); a[0]=(yyvsp[-1].dval);
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
#line 1910 "regparser.c" /* yacc.c:1646  */
    break;

  case 39:
#line 308 "regparser.y" /* yacc.c:1646  */
    {
       double x[1]; double y[1]; double r[2]; double a[1];
       x[0]=(yyvsp[-9].dval); y[0]=(yyvsp[-7].dval); r[0]=(yyvsp[-5].dval); r[1]=(yyvsp[-3].dval); a[0]=(yyvsp[-1].dval);
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
#line 1929 "regparser.c" /* yacc.c:1646  */
    break;

  case 40:
#line 323 "regparser.y" /* yacc.c:1646  */
    {
       double x[1]; double y[1]; double r[2]; double a[1];
       x[0]=(yyvsp[-9].dval); y[0]=(yyvsp[-7].dval); r[0]=(yyvsp[-5].dval); r[1]=(yyvsp[-3].dval); a[0]=(yyvsp[-1].dval);
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
#line 1948 "regparser.c" /* yacc.c:1646  */
    break;

  case 41:
#line 338 "regparser.y" /* yacc.c:1646  */
    {
       double x[1]; double y[1]; double r[2]; double a[1];
       x[0]=(yyvsp[-9].dval); y[0]=(yyvsp[-7].dval); r[0]=(yyvsp[-5].dval); r[1]=(yyvsp[-3].dval); a[0]=(yyvsp[-1].dval);
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
#line 1967 "regparser.c" /* yacc.c:1646  */
    break;

  case 42:
#line 353 "regparser.y" /* yacc.c:1646  */
    {
       double x[1]; double y[1]; double a[2]; double r[2];
       x[0]=(yyvsp[-11].dval); y[0]=(yyvsp[-9].dval); r[0]=(yyvsp[-7].dval); r[1]=(yyvsp[-5].dval);  a[0]=(yyvsp[-3].dval); a[1]=(yyvsp[-1].dval);
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
#line 1990 "regparser.c" /* yacc.c:1646  */
    break;

  case 43:
#line 372 "regparser.y" /* yacc.c:1646  */
    {
       double x[1]; double y[1]; double a[2]; double r[2];
       x[0]=(yyvsp[-11].dval); y[0]=(yyvsp[-9].dval); r[0] = (yyvsp[-7].dval); r[1] = (yyvsp[-5].dval); a[0]=(yyvsp[-3].dval); a[1]=(yyvsp[-1].dval);
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
#line 2013 "regparser.c" /* yacc.c:1646  */
    break;

  case 44:
#line 391 "regparser.y" /* yacc.c:1646  */
    {
       double x[1]; double y[1]; double a[2];
       x[0]=(yyvsp[-7].dval); y[0]=(yyvsp[-5].dval); a[0]=(yyvsp[-3].dval); a[1]=(yyvsp[-1].dval);
       (yyval.my_shape) = regCreateNewWorldShape( regSECTOR, regInclude, x, y, 1, NULL, a, world_coord, 0 );
       if ( (yyval.my_shape) == NULL ) {
 	 my_Gregion = NULL;
	 YYERROR;
       }
     }
#line 2027 "regparser.c" /* yacc.c:1646  */
    break;

  case 45:
#line 401 "regparser.y" /* yacc.c:1646  */
    {
       double x[1]; double y[1]; double a[2];
       x[0]=(yyvsp[-7].dval); y[0]=(yyvsp[-5].dval); a[0]=(yyvsp[-3].dval); a[1]=(yyvsp[-1].dval);
       (yyval.my_shape) = regCreateNewWorldShape( regSECTOR, regExclude, x, y, 1, NULL, a, world_coord, 0 );
       if ( (yyval.my_shape) == NULL ) {
 	 my_Gregion = NULL;
	 YYERROR;
       }
     }
#line 2041 "regparser.c" /* yacc.c:1646  */
    break;

  case 46:
#line 411 "regparser.y" /* yacc.c:1646  */
    {
       double x[1]; double y[1];
       x[0]=(yyvsp[-3].dval); y[0]=(yyvsp[-1].dval);
       (yyval.my_shape) = regCreateNewWorldShape( regPOINT, regInclude, x, y, 1, NULL, NULL, world_coord, 0 );
       if ( (yyval.my_shape) == NULL ) {
 	 my_Gregion = NULL;
	 YYERROR;
       }
     }
#line 2055 "regparser.c" /* yacc.c:1646  */
    break;

  case 47:
#line 421 "regparser.y" /* yacc.c:1646  */
    {
       double x[1]; double y[1];
       x[0]=(yyvsp[-3].dval); y[0]=(yyvsp[-1].dval);
       (yyval.my_shape) = regCreateNewWorldShape( regPOINT, regExclude, x, y, 1, NULL, NULL, world_coord, 0 );
       if ( (yyval.my_shape) == NULL ) {
 	 my_Gregion = NULL;
	 YYERROR;
       }
     }
#line 2069 "regparser.c" /* yacc.c:1646  */
    break;

  case 48:
#line 431 "regparser.y" /* yacc.c:1646  */
    {
       double x[2]; double y[2];
       x[0]=(yyvsp[-7].dval); x[1]=(yyvsp[-3].dval); y[0]=(yyvsp[-5].dval);y[1]=(yyvsp[-1].dval);
       (yyval.my_shape) = regCreateNewWorldShape( regRECTANGLE, regInclude, x, y, 2, NULL, NULL, world_coord, 0 );
       if ( (yyval.my_shape) == NULL ) {
 	 my_Gregion = NULL;
	 YYERROR;
       }
     }
#line 2083 "regparser.c" /* yacc.c:1646  */
    break;

  case 49:
#line 441 "regparser.y" /* yacc.c:1646  */
    {
       double x[2]; double y[2];
       x[0]=(yyvsp[-7].dval); x[1]=(yyvsp[-3].dval); y[0]=(yyvsp[-5].dval); y[1]=(yyvsp[-1].dval);
       (yyval.my_shape) = regCreateNewWorldShape( regRECTANGLE, regExclude, x, y, 2, NULL, NULL, world_coord, 0 );
       if ( (yyval.my_shape) == NULL ) {
 	 my_Gregion = NULL;
	 YYERROR;
       }
     }
#line 2097 "regparser.c" /* yacc.c:1646  */
    break;

  case 50:
#line 451 "regparser.y" /* yacc.c:1646  */
    {
       /* RegLine doesn't work correctly; need to calculate angle */
       double x[2]; double y[2];
       x[0]=(yyvsp[-7].dval); x[1]=(yyvsp[-3].dval); y[0]=(yyvsp[-5].dval);y[1]=(yyvsp[-1].dval);
       (yyval.my_shape) = regCreateNewWorldShape( regRECTANGLE, regInclude, x, y, 2, NULL, NULL, world_coord, 0 );
       if ( (yyval.my_shape) == NULL ) {
 	 my_Gregion = NULL;
	 YYERROR;
       }
     }
#line 2112 "regparser.c" /* yacc.c:1646  */
    break;

  case 51:
#line 462 "regparser.y" /* yacc.c:1646  */
    {
       double x[2]; double y[2];
       x[0]=(yyvsp[-7].dval); x[1]=(yyvsp[-3].dval); y[0]=(yyvsp[-5].dval); y[1]=(yyvsp[-1].dval);
       (yyval.my_shape) = regCreateNewWorldShape( regRECTANGLE, regExclude, x, y, 2, NULL, NULL, world_coord, 0 );
       if ( (yyval.my_shape) == NULL ) {
 	 my_Gregion = NULL;
	 YYERROR;
       }
     }
#line 2126 "regparser.c" /* yacc.c:1646  */
    break;

  case 52:
#line 472 "regparser.y" /* yacc.c:1646  */
    {
       double x[1]; double y[1]; double r[2]; double a[1];
       x[0]=(yyvsp[-9].dval); y[0]=(yyvsp[-7].dval); r[0]=(yyvsp[-5].dval); r[1]=(yyvsp[-3].dval); a[0]=(yyvsp[-1].dval);
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
#line 2145 "regparser.c" /* yacc.c:1646  */
    break;

  case 53:
#line 487 "regparser.y" /* yacc.c:1646  */
    {
       double x[1]; double y[1]; double r[2]; double a[1];
       x[0]=(yyvsp[-9].dval); y[0]=(yyvsp[-7].dval); r[0]=(yyvsp[-5].dval); r[1]=(yyvsp[-3].dval); a[0]=(yyvsp[-1].dval);
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
#line 2164 "regparser.c" /* yacc.c:1646  */
    break;

  case 54:
#line 503 "regparser.y" /* yacc.c:1646  */
    { 
       (yyval.my_shape) = regCreateNewWorldShape( regPOLYGON, regInclude, (yyvsp[-1].PolySide).polyX, (yyvsp[-1].PolySide).polyY, 
			    (yyvsp[-1].PolySide).polyS, NULL, NULL, world_coord, 0 );
       if ( (yyval.my_shape) == NULL ) {
	 my_Gregion = NULL;
	 YYERROR;
       }

       free( (yyvsp[-1].PolySide).polyX );
       free( (yyvsp[-1].PolySide).polyY );
     }
#line 2180 "regparser.c" /* yacc.c:1646  */
    break;

  case 55:
#line 515 "regparser.y" /* yacc.c:1646  */
    { 
       (yyval.my_shape) = regCreateNewWorldShape( regPOLYGON, regExclude, (yyvsp[-1].PolySide).polyX, (yyvsp[-1].PolySide).polyY, 
			    (yyvsp[-1].PolySide).polyS, NULL, NULL, world_coord, 0 );
       free( (yyvsp[-1].PolySide).polyX );
       free( (yyvsp[-1].PolySide).polyY );
       if ( (yyval.my_shape) == NULL ) {
	 my_Gregion = NULL;
	 YYERROR;
       }
     }
#line 2195 "regparser.c" /* yacc.c:1646  */
    break;

  case 56:
#line 529 "regparser.y" /* yacc.c:1646  */
    {
        double x[1]; double y[1]; double r1[2]; double r2[2];
        x[0]=(yyvsp[-11].dval); y[0]=(yyvsp[-9].dval); r1[0]=(yyvsp[-7].dval); r1[1]=(yyvsp[-5].dval); r2[0] = (yyvsp[-3].dval); r2[1] = (yyvsp[-1].dval);
 	    
        fprintf(stderr, "ERROR: Elliptannuli are not yet supported.\n");
        my_Gregion = NULL;
	    YYERROR;
     }
#line 2208 "regparser.c" /* yacc.c:1646  */
    break;

  case 57:
#line 538 "regparser.y" /* yacc.c:1646  */
    {
        double x[1]; double y[1]; double r1[2]; double r2[2];  double a[1];
        x[0]=(yyvsp[-13].dval); y[0]=(yyvsp[-11].dval); r1[0]=(yyvsp[-9].dval); r1[1]=(yyvsp[-7].dval); r2[0] = (yyvsp[-5].dval); r2[1] = (yyvsp[-3].dval); a[0]=(yyvsp[-1].dval);
 	    
        fprintf(stderr, "ERROR: Elliptannuli are not yet supported.\n");
        my_Gregion = NULL;
	    YYERROR;
     }
#line 2221 "regparser.c" /* yacc.c:1646  */
    break;

  case 58:
#line 547 "regparser.y" /* yacc.c:1646  */
    {
        double x[1]; double y[1]; double r1[2]; double r2[2];
        x[0]=(yyvsp[-11].dval); y[0]=(yyvsp[-9].dval); r1[0]=(yyvsp[-7].dval); r1[1]=(yyvsp[-5].dval); r2[0] = (yyvsp[-3].dval); r2[1] = (yyvsp[-1].dval);
 	    
        fprintf(stderr, "ERROR: Elliptannuli are not yet supported.\n");
        my_Gregion = NULL;
	    YYERROR;
     }
#line 2234 "regparser.c" /* yacc.c:1646  */
    break;

  case 59:
#line 556 "regparser.y" /* yacc.c:1646  */
    {
        double x[1]; double y[1]; double r1[2]; double r2[2];  double a[1];
        x[0]=(yyvsp[-13].dval); y[0]=(yyvsp[-11].dval); r1[0]=(yyvsp[-9].dval); r1[1]=(yyvsp[-7].dval); r2[0] = (yyvsp[-5].dval); r2[1] = (yyvsp[-3].dval); a[0]=(yyvsp[-1].dval);
 	    
        fprintf(stderr, "ERROR: Elliptannuli are not yet supported.\n");
        my_Gregion = NULL;
	    YYERROR;
     }
#line 2247 "regparser.c" /* yacc.c:1646  */
    break;

  case 60:
#line 568 "regparser.y" /* yacc.c:1646  */
    {
        double x[1]; double y[1]; double r[2]; 
        x[0]=(yyvsp[-7].dval); y[0]=(yyvsp[-5].dval); r[0]=(yyvsp[-3].dval); r[1]=(yyvsp[-1].dval); 
 	    
        fprintf(stderr, "ERROR: Diamonds are not yet supported.\n");
        my_Gregion = NULL;
	    YYERROR;
     }
#line 2260 "regparser.c" /* yacc.c:1646  */
    break;

  case 61:
#line 577 "regparser.y" /* yacc.c:1646  */
    {
        double x[1]; double y[1]; double r[2]; 
        x[0]=(yyvsp[-7].dval); y[0]=(yyvsp[-5].dval); r[0]=(yyvsp[-3].dval); r[1]=(yyvsp[-1].dval); 
 	    
        fprintf(stderr, "ERROR: Diamonds are not yet supported.\n");
        my_Gregion = NULL;
	    YYERROR;
     }
#line 2273 "regparser.c" /* yacc.c:1646  */
    break;

  case 62:
#line 586 "regparser.y" /* yacc.c:1646  */
    {
        double x[1]; double y[1]; double r[2]; double a[1];
        x[0]=(yyvsp[-9].dval); y[0]=(yyvsp[-7].dval); r[0]=(yyvsp[-5].dval); r[1]=(yyvsp[-3].dval); a[0] = (yyvsp[-1].dval);
 	    
        fprintf(stderr, "ERROR: Diamonds are not yet supported.\n");
        my_Gregion = NULL;
	    YYERROR;
     }
#line 2286 "regparser.c" /* yacc.c:1646  */
    break;

  case 63:
#line 595 "regparser.y" /* yacc.c:1646  */
    {
        double x[1]; double y[1]; double r[2]; double a[1];
        x[0]=(yyvsp[-9].dval); y[0]=(yyvsp[-7].dval); r[0]=(yyvsp[-5].dval); r[1]=(yyvsp[-3].dval); a[0] = (yyvsp[-1].dval);
 	    
        fprintf(stderr, "ERROR: Diamonds are not yet supported.\n");
        my_Gregion = NULL;
	    YYERROR;
     }
#line 2299 "regparser.c" /* yacc.c:1646  */
    break;

  case 64:
#line 607 "regparser.y" /* yacc.c:1646  */
    {
     (yyval.PolySide) = (yyvsp[-4].PolySide);
     (yyval.PolySide).polyS += 1;
     (yyval.PolySide).polyX = (double*)realloc( (yyval.PolySide).polyX, (yyval.PolySide).polyS *sizeof(double));
     (yyval.PolySide).polyY = (double*)realloc( (yyval.PolySide).polyY, (yyval.PolySide).polyS *sizeof(double));
     (yyval.PolySide).polyX[(yyval.PolySide).polyS -1] = (yyvsp[-2].dval);
     (yyval.PolySide).polyY[(yyval.PolySide).polyS -1] = (yyvsp[0].dval);
     }
#line 2312 "regparser.c" /* yacc.c:1646  */
    break;

  case 65:
#line 615 "regparser.y" /* yacc.c:1646  */
    {
     (yyval.PolySide).polyS = 1;
     (yyval.PolySide).polyX = (double*)calloc(1,sizeof(double));
     (yyval.PolySide).polyY = (double*)calloc(1,sizeof(double));
     (yyval.PolySide).polyX[0] = (yyvsp[-2].dval);
     (yyval.PolySide).polyY[0] = (yyvsp[0].dval);
   }
#line 2324 "regparser.c" /* yacc.c:1646  */
    break;


#line 2328 "regparser.c" /* yacc.c:1646  */
      default: break;
    }
  /* User semantic actions sometimes alter yychar, and that requires
     that yytoken be updated with the new translation.  We take the
     approach of translating immediately before every use of yytoken.
     One alternative is translating here after every semantic action,
     but that translation would be missed if the semantic action invokes
     YYABORT, YYACCEPT, or YYERROR immediately after altering yychar or
     if it invokes YYBACKUP.  In the case of YYABORT or YYACCEPT, an
     incorrect destructor might then be invoked immediately.  In the
     case of YYERROR or YYBACKUP, subsequent parser actions might lead
     to an incorrect destructor call or verbose syntax error message
     before the lookahead is translated.  */
  YY_SYMBOL_PRINT ("-> $$ =", yyr1[yyn], &yyval, &yyloc);

  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);

  *++yyvsp = yyval;

  /* Now 'shift' the result of the reduction.  Determine what state
     that goes to, based on the state we popped back to and the rule
     number reduced by.  */

  yyn = yyr1[yyn];

  yystate = yypgoto[yyn - YYNTOKENS] + *yyssp;
  if (0 <= yystate && yystate <= YYLAST && yycheck[yystate] == *yyssp)
    yystate = yytable[yystate];
  else
    yystate = yydefgoto[yyn - YYNTOKENS];

  goto yynewstate;


/*--------------------------------------.
| yyerrlab -- here on detecting error.  |
`--------------------------------------*/
yyerrlab:
  /* Make sure we have latest lookahead translation.  See comments at
     user semantic actions for why this is necessary.  */
  yytoken = yychar == YYEMPTY ? YYEMPTY : YYTRANSLATE (yychar);

  /* If not already recovering from an error, report this error.  */
  if (!yyerrstatus)
    {
      ++yynerrs;
#if ! YYERROR_VERBOSE
      yyerror (YY_("syntax error"));
#else
# define YYSYNTAX_ERROR yysyntax_error (&yymsg_alloc, &yymsg, \
                                        yyssp, yytoken)
      {
        char const *yymsgp = YY_("syntax error");
        int yysyntax_error_status;
        yysyntax_error_status = YYSYNTAX_ERROR;
        if (yysyntax_error_status == 0)
          yymsgp = yymsg;
        else if (yysyntax_error_status == 1)
          {
            if (yymsg != yymsgbuf)
              YYSTACK_FREE (yymsg);
            yymsg = (char *) YYSTACK_ALLOC (yymsg_alloc);
            if (!yymsg)
              {
                yymsg = yymsgbuf;
                yymsg_alloc = sizeof yymsgbuf;
                yysyntax_error_status = 2;
              }
            else
              {
                yysyntax_error_status = YYSYNTAX_ERROR;
                yymsgp = yymsg;
              }
          }
        yyerror (yymsgp);
        if (yysyntax_error_status == 2)
          goto yyexhaustedlab;
      }
# undef YYSYNTAX_ERROR
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

  /* Do not reclaim the symbols of the rule whose action triggered
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
  yyerrstatus = 3;      /* Each real token shifted decrements this.  */

  for (;;)
    {
      yyn = yypact[yystate];
      if (!yypact_value_is_default (yyn))
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

  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  *++yyvsp = yylval;
  YY_IGNORE_MAYBE_UNINITIALIZED_END


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

#if !defined yyoverflow || YYERROR_VERBOSE
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
    {
      /* Make sure we have latest lookahead translation.  See comments at
         user semantic actions for why this is necessary.  */
      yytoken = YYTRANSLATE (yychar);
      yydestruct ("Cleanup: discarding lookahead",
                  yytoken, &yylval);
    }
  /* Do not reclaim the symbols of the rule whose action triggered
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
  return yyresult;
}
#line 623 "regparser.y" /* yacc.c:1906  */



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

