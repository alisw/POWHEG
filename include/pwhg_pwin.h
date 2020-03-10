c -*- Fortran -*-

      integer maxnum,maxstrings,maxlin,maxkey
      parameter (maxnum=250,maxstrings=20,maxlin=100,maxkey=20)
      integer pwin_numvalues,pwin_numstrings
      character * (maxkey) pwin_keywords(maxnum)
      real * 8 pwin_values(maxnum)
      logical pwin_used(maxnum)
      integer pwin_stringptr(maxnum)
      character * (maxlin) pwin_strings(maxstrings)
      common/pwhg_pwin/pwin_keywords,pwin_stringptr,pwin_values,
     1     pwin_strings,pwin_used,pwin_numvalues,pwin_numstrings
