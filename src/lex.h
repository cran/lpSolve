#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
# define U(x) x
# define NLSTATE yyprevious=YYNEWLINE
# define BEGIN yybgin = yysvec + 1 +
# define INITIAL 0
# define YYLERR yysvec
# define YYSTATE (yyestate-yysvec-1)
# define YYOPTIM 1
# ifndef YYLMAX 
# define YYLMAX BUFSIZ
# endif 
#ifndef __cplusplus
# define output(c) (void)putc(c,yyout)
#else
# define lex_output(c) (void)putc(c,yyout)
#endif

#if defined(__cplusplus) || defined(__STDC__)

#if defined(__cplusplus) && defined(__EXTERN_C__)
extern "C" {
#endif
	int yyback(int *, int);
	int yyinput(void);
	int yylook(void);
	void yyoutput(int);
	int yyracc(int);
	int yyreject(void);
	void yyunput(int);
	int yylex(void);
#ifdef YYLEX_E
	void yywoutput(wchar_t);
	wchar_t yywinput(void);
#endif
#ifndef yyless
	int yyless(int);
#endif
#ifndef yywrap
	int yywrap(void);
#endif
#ifdef LEXDEBUG
	void allprint(char);
	void sprint(char *);
#endif
#if defined(__cplusplus) && defined(__EXTERN_C__)
}
#endif

#ifdef __cplusplus
extern "C" {
#endif
	void exit(int);
#ifdef __cplusplus
}
#endif

#endif
# define unput(c) {yytchar= (c);if(yytchar=='\n')yylineno--;*yysptr++=yytchar;}
# define yymore() (yymorfg=1)
#ifndef __cplusplus
# define input() (((yytchar=yysptr>yysbuf?U(*--yysptr):getc(yyin))==10?(yylineno++,yytchar):yytchar)==EOF?0:yytchar)
#else
# define lex_input() (((yytchar=yysptr>yysbuf?U(*--yysptr):getc(yyin))==10?(yylineno++,yytchar):yytchar)==EOF?0:yytchar)
#endif
#define ECHO fprintf(yyout, "%s",yytext)
# define REJECT { nstr = yyreject(); goto yyfussy;}
int yyleng;
#define YYISARRAY
char yytext[YYLMAX];
int yymorfg;
extern char *yysptr, yysbuf[];
int yytchar;
FILE *yyin, *yyout;
extern int yylineno;
struct yysvf { 
	struct yywork *yystoff;
	struct yysvf *yyother;
	int *yystops;};
struct yysvf *yyestate;
extern struct yysvf yysvec[], *yybgin;
# define COMMENT 2
# define YYNEWLINE 10
int yylex(){
int nstr;
/* int nstr; extern int yyprevious; */
#ifdef __cplusplus
/* to avoid CC and lint complaining yyfussy not being used ...*/
static int __lex_hack = 0;
if (__lex_hack) goto yyfussy;
#endif
yyin = stdin; yyout = stdout;
while((nstr = yylook()) >= 0)
/* yyfussy: switch(nstr){ */
switch(nstr){
case 0:
if(yywrap()) return(0); break;
case 1:

# line 13 "lex.l.old"
{
  BEGIN COMMENT;
}
break;
case 2:

# line 17 "lex.l.old"
{
  BEGIN INITIAL;
}
break;
case 3:

# line 21 "lex.l.old"
{
}
break;
case 4:

# line 24 "lex.l.old"
{
}
break;
case 5:

# line 27 "lex.l.old"
{
}
break;
case 6:

# line 30 "lex.l.old"
{
  return(COMMA);
}
break;
case 7:

# line 34 "lex.l.old"
{
  return(MINIMISE);
}
break;
case 8:

# line 38 "lex.l.old"
{
  return(MAXIMISE);
}
break;
case 9:

# line 42 "lex.l.old"
{
  f = atof((char *)yytext);
  return(CONS);
}
break;
case 10:

# line 47 "lex.l.old"
{
  Sign = 0;
  for(x = 0; x < yyleng; x++)
    if(yytext[x] == '-' || yytext[x] == '+')
      Sign = (Sign == (yytext[x] == '+'));
  return (SIGN);
  /* Sign is TRUE if the sign-string
     represents a '-'. Otherwise Sign
     is FALSE */
}
break;
case 11:

# line 58 "lex.l.old"
{
  Within_int_decl = TRUE;
  return(VAR);
}
break;
case 12:

# line 63 "lex.l.old"
{
  strcpy(Last_var, (char *)yytext);
  return(VAR);
}
break;
case 13:

# line 68 "lex.l.old"
{
  return (COLON);
}
break;
case 14:

# line 72 "lex.l.old"
{
  return(AR_M_OP);
}
break;
case 15:

# line 76 "lex.l.old"
{
  return(RE_OP);
}
break;
case 16:

# line 80 "lex.l.old"
{
  Within_int_decl = FALSE;
  return(END_C);
}
break;
case 17:

# line 85 "lex.l.old"
{
  report(NULL, CRITICALSTOP,"LEX ERROR : %s lineno %d \n" ,yytext,yylineno);
}
break;
case -1:
break;
default:
(void)fprintf(yyout,"bad switch yylook %d",nstr);
} return(0); }
/* end of yylex */
int yyvstop[] = {
0,

15,
0, 

15,
0, 

15,
0, 

15,
0, 

17,
0, 

5,
10,
17,
0, 

5,
10,
0, 

14,
17,
0, 

10,
17,
0, 

6,
17,
0, 

17,
0, 

17,
0, 

9,
17,
0, 

13,
17,
0, 

16,
17,
0, 

15,
17,
0, 

15,
17,
0, 

12,
17,
0, 

12,
17,
0, 

12,
17,
0, 

3,
0, 

4,
0, 

3,
0, 

10,
0, 

9,
0, 

1,
0, 

9,
0, 

15,
0, 

12,
0, 

12,
0, 

12,
0, 

12,
0, 

2,
0, 

9,
0, 

11,
12,
0, 

12,
0, 

12,
0, 

12,
0, 

8,
0, 

7,
0, 
0};
# define YYTYPE unsigned char
struct yywork { YYTYPE verify, advance; } yycrank[] = {
{0,0},	{0,0},	{1,5},	{0,0},	
{0,0},	{0,0},	{0,0},	{0,0},	
{0,0},	{0,0},	{1,6},	{1,7},	
{9,24},	{9,24},	{0,0},	{0,0},	
{0,0},	{6,7},	{6,7},	{0,0},	
{0,0},	{0,0},	{0,0},	{0,0},	
{0,0},	{0,0},	{0,0},	{0,0},	
{0,0},	{0,0},	{0,0},	{0,0},	
{24,24},	{24,24},	{0,0},	{9,24},	
{1,5},	{0,0},	{0,0},	{0,0},	
{6,7},	{0,0},	{0,0},	{1,8},	
{1,9},	{1,10},	{4,23},	{1,11},	
{1,12},	{1,13},	{12,26},	{6,24},	
{2,8},	{6,24},	{2,10},	{24,24},	
{2,11},	{2,12},	{23,36},	{1,14},	
{1,15},	{1,16},	{1,17},	{16,30},	
{3,21},	{0,0},	{1,18},	{1,18},	
{2,14},	{2,15},	{1,18},	{2,17},	
{3,21},	{3,22},	{1,19},	{39,43},	
{41,44},	{42,45},	{1,20},	{1,18},	
{0,0},	{0,0},	{0,0},	{0,0},	
{0,0},	{1,18},	{0,0},	{0,0},	
{0,0},	{1,18},	{19,33},	{0,0},	
{25,29},	{20,34},	{35,42},	{33,40},	
{0,0},	{0,0},	{3,21},	{34,41},	
{0,0},	{20,35},	{0,0},	{0,0},	
{0,0},	{3,23},	{3,21},	{0,0},	
{0,0},	{0,0},	{0,0},	{3,21},	
{11,25},	{11,25},	{11,25},	{11,25},	
{11,25},	{11,25},	{11,25},	{11,25},	
{11,25},	{11,25},	{19,33},	{3,21},	
{25,29},	{20,34},	{35,42},	{33,40},	
{3,21},	{3,21},	{0,0},	{34,41},	
{3,21},	{20,35},	{0,0},	{0,0},	
{3,21},	{0,0},	{0,0},	{0,0},	
{3,21},	{3,21},	{0,0},	{0,0},	
{0,0},	{0,0},	{0,0},	{3,21},	
{0,0},	{0,0},	{13,27},	{3,21},	
{13,28},	{13,28},	{13,28},	{13,28},	
{13,28},	{13,28},	{13,28},	{13,28},	
{13,28},	{13,28},	{37,38},	{37,38},	
{37,38},	{37,38},	{37,38},	{37,38},	
{37,38},	{37,38},	{37,38},	{37,38},	
{0,0},	{13,29},	{0,0},	{0,0},	
{0,0},	{18,31},	{18,31},	{18,31},	
{18,31},	{18,31},	{0,0},	{0,0},	
{0,0},	{0,0},	{0,0},	{0,0},	
{18,31},	{18,31},	{18,31},	{18,31},	
{18,31},	{18,31},	{18,31},	{18,31},	
{18,31},	{18,31},	{18,31},	{18,31},	
{0,0},	{0,0},	{18,32},	{0,0},	
{0,0},	{13,29},	{18,31},	{18,31},	
{18,31},	{18,31},	{18,31},	{18,31},	
{18,31},	{18,31},	{18,31},	{18,31},	
{18,31},	{18,31},	{18,31},	{18,31},	
{18,31},	{18,31},	{18,31},	{18,31},	
{18,31},	{18,31},	{18,31},	{18,31},	
{18,31},	{18,31},	{18,31},	{18,31},	
{18,31},	{18,31},	{0,0},	{18,31},	
{18,31},	{18,31},	{0,0},	{18,31},	
{18,31},	{18,31},	{18,31},	{18,31},	
{18,31},	{18,31},	{18,31},	{18,31},	
{18,31},	{18,31},	{18,31},	{18,31},	
{18,31},	{18,31},	{18,31},	{18,31},	
{18,31},	{18,31},	{18,31},	{18,31},	
{18,31},	{18,31},	{18,31},	{18,31},	
{18,31},	{18,31},	{0,0},	{18,31},	
{18,31},	{29,37},	{0,0},	{29,37},	
{0,0},	{0,0},	{29,38},	{29,38},	
{29,38},	{29,38},	{29,38},	{29,38},	
{29,38},	{29,38},	{29,38},	{29,38},	
{32,39},	{32,39},	{32,39},	{32,39},	
{32,39},	{0,0},	{0,0},	{0,0},	
{0,0},	{0,0},	{0,0},	{32,39},	
{32,39},	{32,39},	{32,39},	{32,39},	
{32,39},	{32,39},	{32,39},	{32,39},	
{32,39},	{32,39},	{32,39},	{0,0},	
{0,0},	{0,0},	{0,0},	{0,0},	
{0,0},	{32,39},	{32,39},	{32,39},	
{32,39},	{32,39},	{32,39},	{32,39},	
{32,39},	{32,39},	{32,39},	{32,39},	
{32,39},	{32,39},	{32,39},	{32,39},	
{32,39},	{32,39},	{32,39},	{32,39},	
{32,39},	{32,39},	{32,39},	{32,39},	
{32,39},	{32,39},	{32,39},	{32,39},	
{32,39},	{0,0},	{32,39},	{32,39},	
{32,39},	{0,0},	{32,39},	{32,39},	
{32,39},	{32,39},	{32,39},	{32,39},	
{32,39},	{32,39},	{32,39},	{32,39},	
{32,39},	{32,39},	{32,39},	{32,39},	
{32,39},	{32,39},	{32,39},	{32,39},	
{32,39},	{32,39},	{32,39},	{32,39},	
{32,39},	{32,39},	{32,39},	{32,39},	
{32,39},	{0,0},	{32,39},	{32,39},	
{0,0}};
struct yysvf yysvec[] = {
{ 0,	0,	0},
{ yycrank+-1,	0,		yyvstop+1},
{ yycrank+-10,	yysvec+1,	yyvstop+3},
{ yycrank+-63,	0,		yyvstop+5},
{ yycrank+-4,	yysvec+3,	yyvstop+7},
{ yycrank+0,	0,		yyvstop+9},
{ yycrank+8,	0,		yyvstop+11},
{ yycrank+0,	yysvec+6,	yyvstop+15},
{ yycrank+0,	0,		yyvstop+18},
{ yycrank+3,	yysvec+6,	yyvstop+21},
{ yycrank+0,	0,		yyvstop+24},
{ yycrank+64,	0,		yyvstop+27},
{ yycrank+8,	0,		yyvstop+29},
{ yycrank+104,	0,		yyvstop+31},
{ yycrank+0,	0,		yyvstop+34},
{ yycrank+0,	0,		yyvstop+37},
{ yycrank+2,	0,		yyvstop+40},
{ yycrank+0,	0,		yyvstop+43},
{ yycrank+142,	0,		yyvstop+46},
{ yycrank+12,	yysvec+18,	yyvstop+49},
{ yycrank+28,	yysvec+18,	yyvstop+52},
{ yycrank+0,	0,		yyvstop+55},
{ yycrank+0,	0,		yyvstop+57},
{ yycrank+11,	0,		yyvstop+59},
{ yycrank+23,	yysvec+6,	yyvstop+61},
{ yycrank+23,	yysvec+11,	yyvstop+63},
{ yycrank+0,	0,		yyvstop+65},
{ yycrank+0,	yysvec+11,	0},	
{ yycrank+0,	yysvec+13,	yyvstop+67},
{ yycrank+226,	0,		0},	
{ yycrank+0,	0,		yyvstop+69},
{ yycrank+0,	yysvec+18,	yyvstop+71},
{ yycrank+249,	0,		0},	
{ yycrank+11,	yysvec+18,	yyvstop+73},
{ yycrank+11,	yysvec+18,	yyvstop+75},
{ yycrank+16,	yysvec+18,	yyvstop+77},
{ yycrank+0,	0,		yyvstop+79},
{ yycrank+114,	0,		0},	
{ yycrank+0,	yysvec+37,	yyvstop+81},
{ yycrank+13,	yysvec+32,	0},	
{ yycrank+0,	yysvec+18,	yyvstop+83},
{ yycrank+18,	yysvec+18,	yyvstop+86},
{ yycrank+19,	yysvec+18,	yyvstop+88},
{ yycrank+0,	0,		yyvstop+90},
{ yycrank+0,	0,		yyvstop+92},
{ yycrank+0,	0,		yyvstop+94},
{ 0,	0,	0}};
struct yywork *yytop = yycrank+375;
struct yysvf *yybgin = yysvec+1;
char yymatch[] = {
  0,   1,   1,   1,   1,   1,   1,   1, 
  1,   9,  10,   1,   1,   1,   1,   1, 
  1,   1,   1,   1,   1,   1,   1,   1, 
  1,   1,   1,   1,   1,   1,   1,   1, 
  9,   1,   1,  35,  35,  35,  35,  35, 
  1,   1,   1,  43,   1,  43,  35,  35, 
 48,  48,  48,  48,  48,  48,  48,  48, 
 48,  48,   1,   1,  60,   1,  60,   1, 
 35,  65,  66,  66,  66,  69,  66,  66, 
 66,  73,  66,  66,  66,  77,  78,  66, 
 66,  66,  66,  66,  84,  66,  66,  66, 
 88,  66,  66,  35,   1,  35,  35,  35, 
  1,  65,  66,  66,  66,  69,  66,  66, 
 66,  73,  66,  66,  66,  77,  78,  66, 
 66,  66,  66,  66,  84,  66,  66,  66, 
 88,  66,  66,  35,   1,  35,  35,   1, 
  1,   1,   1,   1,   1,   1,   1,   1, 
  1,   1,   1,   1,   1,   1,   1,   1, 
  1,   1,   1,   1,   1,   1,   1,   1, 
  1,   1,   1,   1,   1,   1,   1,   1, 
  1,   1,   1,   1,   1,   1,   1,   1, 
  1,   1,   1,   1,   1,   1,   1,   1, 
  1,   1,   1,   1,   1,   1,   1,   1, 
  1,   1,   1,   1,   1,   1,   1,   1, 
  1,   1,   1,   1,   1,   1,   1,   1, 
  1,   1,   1,   1,   1,   1,   1,   1, 
  1,   1,   1,   1,   1,   1,   1,   1, 
  1,   1,   1,   1,   1,   1,   1,   1, 
  1,   1,   1,   1,   1,   1,   1,   1, 
  1,   1,   1,   1,   1,   1,   1,   1, 
  1,   1,   1,   1,   1,   1,   1,   1, 
  1,   1,   1,   1,   1,   1,   1,   1, 
0};
char yyextra[] = {
0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,
0};
/*	Copyright (c) 1989 AT&T	*/
/*	  All Rights Reserved  	*/

/*	THIS IS UNPUBLISHED PROPRIETARY SOURCE CODE OF AT&T	*/
/*	The copyright notice above does not evidence any   	*/
/*	actual or intended publication of such source code.	*/

/* #pragma ident	"@(#)ncform	6.12	97/12/08 SMI" */

int yylineno =1;
# define YYU(x) x
# define NLSTATE yyprevious=YYNEWLINE
struct yysvf *yylstate [YYLMAX], **yylsp, **yyolsp;
char yysbuf[YYLMAX];
char *yysptr = yysbuf;
int *yyfnd;
extern struct yysvf *yyestate;
int yyprevious = YYNEWLINE;
#if defined(__cplusplus) || defined(__STDC__)
int yylook(void)
#else
yylook()
#endif
{
	register struct yysvf *yystate, **lsp;
	register struct yywork *yyt;
	struct yysvf *yyz;
	int yych, yyfirst;
	struct yywork *yyr;
# ifdef LEXDEBUG
	int debug;
# endif
	char *yylastch;
	/* start off machines */
# ifdef LEXDEBUG
	debug = 0;
# endif
 yyin = stdin; yyout = stdout;
	yyfirst=1;
	if (!yymorfg)
		yylastch = yytext;
	else {
		yymorfg=0;
		yylastch = yytext+yyleng;
		}
	for(;;){
		lsp = yylstate;
		yyestate = yystate = yybgin;
		if (yyprevious==YYNEWLINE) yystate++;
		for (;;){
# ifdef LEXDEBUG
			if(debug)fprintf(yyout,"state %d\n",yystate-yysvec-1);
# endif
			yyt = yystate->yystoff;
			if(yyt == yycrank && !yyfirst){  /* may not be any transitions */
				yyz = yystate->yyother;
				if(yyz == 0)break;
				if(yyz->yystoff == yycrank)break;
				}
#ifndef __cplusplus
			*yylastch++ = yych = input();
#else
			*yylastch++ = yych = lex_input();
#endif
#ifdef YYISARRAY
			if(yylastch > &yytext[YYLMAX]) {
				fprintf(yyout,"Input string too long, limit %d\n",YYLMAX);
				exit(1);
			}
#else
			if (yylastch >= &yytext[ yytextsz ]) {
				int	x = yylastch - yytext;

				yytextsz += YYTEXTSZINC;
				if (yytext == yy_tbuf) {
				    yytext = (char *) malloc(yytextsz);
				    memcpy(yytext, yy_tbuf, sizeof (yy_tbuf));
				}
				else
				    yytext = (char *) realloc(yytext, yytextsz);
				if (!yytext) {
				    fprintf(yyout,
					"Cannot realloc yytext\n");
				    exit(1);
				}
				yylastch = yytext + x;
			}
#endif
			yyfirst=0;
		tryagain:
# ifdef LEXDEBUG
			if(debug){
				fprintf(yyout,"char ");
				allprint(yych);
				putchar('\n');
				}
# endif
			yyr = yyt;
			if ( (uintptr_t)yyt > (uintptr_t)yycrank){
				yyt = yyr + yych;
				if (yyt <= yytop && yyt->verify+yysvec == yystate){
					if(yyt->advance+yysvec == YYLERR)	/* error transitions */
						{unput(*--yylastch);break;}
					*lsp++ = yystate = yyt->advance+yysvec;
					if(lsp > &yylstate[YYLMAX]) {
						fprintf(yyout,"Input string too long, limit %d\n",YYLMAX);
						exit(1);
					}
					goto contin;
					}
				}
# ifdef YYOPTIM
			else if((uintptr_t)yyt < (uintptr_t)yycrank) {	/* r < yycrank */
				yyt = yyr = yycrank+(yycrank-yyt);
# ifdef LEXDEBUG
				if(debug)fprintf(yyout,"compressed state\n");
# endif
				yyt = yyt + yych;
				if(yyt <= yytop && yyt->verify+yysvec == yystate){
					if(yyt->advance+yysvec == YYLERR)	/* error transitions */
						{unput(*--yylastch);break;}
					*lsp++ = yystate = yyt->advance+yysvec;
					if(lsp > &yylstate[YYLMAX]) {
						fprintf(yyout,"Input string too long, limit %d\n",YYLMAX);
						exit(1);
					}
					goto contin;
					}
				yyt = yyr + YYU(yymatch[yych]);
# ifdef LEXDEBUG
				if(debug){
					fprintf(yyout,"try fall back character ");
					allprint(YYU(yymatch[yych]));
					putchar('\n');
					}
# endif
				if(yyt <= yytop && yyt->verify+yysvec == yystate){
					if(yyt->advance+yysvec == YYLERR)	/* error transition */
						{unput(*--yylastch);break;}
					*lsp++ = yystate = yyt->advance+yysvec;
					if(lsp > &yylstate[YYLMAX]) {
						fprintf(yyout,"Input string too long, limit %d\n",YYLMAX);
						exit(1);
					}
					goto contin;
					}
				}
			if ((yystate = yystate->yyother) && (yyt= yystate->yystoff) != yycrank){
# ifdef LEXDEBUG
				if(debug)fprintf(yyout,"fall back to state %d\n",yystate-yysvec-1);
# endif
				goto tryagain;
				}
# endif
			else
				{unput(*--yylastch);break;}
		contin:
# ifdef LEXDEBUG
			if(debug){
				fprintf(yyout,"state %d char ",yystate-yysvec-1);
				allprint(yych);
				putchar('\n');
				}
# endif
			;
			}
# ifdef LEXDEBUG
		if(debug){
			fprintf(yyout,"stopped at %d with ",*(lsp-1)-yysvec-1);
			allprint(yych);
			putchar('\n');
			}
# endif
		while (lsp-- > yylstate){
			*yylastch-- = 0;
			if (*lsp != 0 && (yyfnd= (*lsp)->yystops) && *yyfnd > 0){
				yyolsp = lsp;
				if(yyextra[*yyfnd]){		/* must backup */
					while(yyback((*lsp)->yystops,-*yyfnd) != 1 && lsp > yylstate){
						lsp--;
						unput(*yylastch--);
						}
					}
				yyprevious = YYU(*yylastch);
				yylsp = lsp;
				yyleng = yylastch-yytext+1;
				yytext[yyleng] = 0;
# ifdef LEXDEBUG
				if(debug){
					fprintf(yyout,"\nmatch ");
					sprint(yytext);
					fprintf(yyout," action %d\n",*yyfnd);
					}
# endif
				return(*yyfnd++);
				}
			unput(*yylastch);
			}
		if (yytext[0] == 0  /* && feof(yyin) */)
			{
			yysptr=yysbuf;
			return(0);
			}
#ifndef __cplusplus
		yyprevious = yytext[0] = input();
		if (yyprevious>0)
			output(yyprevious);
#else
		yyprevious = yytext[0] = lex_input();
		if (yyprevious>0)
			lex_output(yyprevious);
#endif
		yylastch=yytext;
# ifdef LEXDEBUG
		if(debug)putchar('\n');
# endif
		}
	}
#if defined(__cplusplus) || defined(__STDC__)
int yyback(int *p, int m)
#else
yyback(p, m)
	int *p;
#endif
{
	if (p==0) return(0);
	while (*p) {
		if (*p++ == m)
			return(1);
	}
	return(0);
}
	/* the following are only used in the lex library */
#if defined(__cplusplus) || defined(__STDC__)
int yyinput(void)
#else
yyinput()
#endif
{
#ifndef __cplusplus
	return(input());
#else
	return(lex_input());
#endif
	}
#if defined(__cplusplus) || defined(__STDC__)
void yyoutput(int c)
#else
yyoutput(c)
  int c; 
#endif
{
yyin = stdin; yyout = stdout;
#ifndef __cplusplus
	output(c);
#else
	lex_output(c);
#endif
	}
#if defined(__cplusplus) || defined(__STDC__)
void yyunput(int c)
#else
yyunput(c)
   int c; 
#endif
{
	unput(c);
	}
