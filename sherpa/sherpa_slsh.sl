% 
%  Copyright (C) 2007-2009  Smithsonian Astrophysical Observatory
%
%
%  This program is free software; you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation; either version 3 of the License, or
%  (at your option) any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License along
%  with this program; if not, write to the Free Software Foundation, Inc.,
%  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
%


_traceback = 0;

private variable _prompt_ctr = 0;
private define _sherpa_prompt(flag) {
  if (flag)
     return "... ";
  _prompt_ctr++;
  return sprintf ("sherpa-%d> ", _prompt_ctr);
}


slsh_set_prompt_hook(&_sherpa_prompt);

require("crates");
require("chips_hlui");

% load resource file- don't display warning/errors in load

set_preference_autoload(1);

% Now load Sherpa, after ChIPS and CRATES
require("sherpa");

% wrap ahelp
define ahelp()
{
   if (_NARGS != 0)
   {
      variable i, buf, cmd;
      variable args = __pop_args (_NARGS);
      buf = "";
      for (i = 0; i < _NARGS; i++)
      {
         buf = sprintf("%s %s", buf, string(args[i].value));
      }
      if (is_substr(buf, "-x ")) 
         cmd = sprintf("ahelp %s", buf);
      else
         cmd = sprintf("ahelp -x '/sl[.]*/' %s", buf);
      () = system(cmd);
   }
   else
   {
      () = system("ahelp");
   }
}


% sherpa script function

#if ( is_defined ("slsh_interactive") )
define script()
{

#if (_slang_version >= 20101)
    variable file = qualifier ("filename","sherpa.log");
    variable clob = qualifier ("clobber", 0); 
#else
    variable file = "sherpa.log", clob = 0;   
#endif

  if (_NARGS > 2) {
       () = fprintf (stderr, "Usage: script( [filename=\"sherpa.log\", clobber=0] )\n");
       return;
  }

  if (_NARGS == 1) {
       file = ();
  }
  else if (_NARGS == 2) {
       (file, clob) = ();
  }
  
  if (stat_file(file) != NULL and clob == 0) {
       () = fprintf (stderr, "script file \"%s\" exists and clobber is not set.\n", file);
       return;
  }

  save_input (file);

}
#endif
