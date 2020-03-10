#!/bin/bash

svnversiondir=`dirname $0`

function svninfo {
if which svn > /dev/null && svn info > /dev/null
then
# only if svn is a command, and current directory is under svn control
url=`svn info | grep -e '^URL'`
revision=`svn info | grep -e '^Revision:'| sed 's/Revision: //'`
else
url=' no svn repository found'
revision=''
fi
}

currdir=`pwd`

> $currdir/svnversion.txt

if [ a"$*" = a ]
then
    dirs="$currdir ../"
else
    dirs="$*"
fi

for dir in $dirs
do

cd $dir

svninfo

pwd >> $currdir/svnversion.txt

echo $url >> $currdir/svnversion.txt

if [ a$revision != a ]
then
echo Rev.$revision >>  $currdir/svnversion.txt

svn status | grep -e '^\?.*\.[fFch]$\|^[MA]'

modlines=`svn status | grep -e '^\?.*\.[fFch]$\|^[MA]'|wc|sed 's/  */ /g'|cut -d' ' -f2`

if [ $modlines = 0 ]
then
    echo "clean version" >> $currdir/svnversion.txt
else
    echo "Warning: not a clean version:"  >> $currdir/svnversion.txt
    svn status | grep -e '^\?.*\.[fFch]$\|^[MA]' | sort >> $currdir/svnversion.txt
fi

fi

cd $currdir

done

gfortran $svnversiondir/svnversion.f -o svnversion

./svnversion

\rm svnversion
\rm svnversion.txt

if ! [ -e svn.version ] || ! cmp svnversion.tmp svn.version > /dev/null 2>&1
then
    \mv svnversion.tmp svn.version
else
\rm svnversion.tmp
fi
