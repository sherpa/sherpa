/*_C_INSERT_SAO_COPYRIGHT_HERE_(2007)_*/
/*_C_INSERT_GPL_LICENSE_HERE_*/

                                                
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










