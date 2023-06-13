#!/usr/bin/sh
export STRVERSION=$1
MAJOR="$(echo $STRVERSION | cut -d'.' -f1)"
MINOR="$(echo $STRVERSION | cut -d'.' -f2)"
MICRO="$(echo $STRVERSION | cut -d'.' -f3)"

error () {
    echo "$@" 1>&2
}

if [ -z $MAJOR ]
then
    error "Major version not specified!"
fi

if [ -z $MINOR ]
then 
    MINOR="0"
fi

if [ -z $MICRO ]
then 
    MICRO="0"
fi

export PYVERSION="\"${MAJOR}\", \"${MINOR}\", \"${MICRO}\""
export PYPROJECTTOMLVERSION="${MAJOR}.${MINOR}.${MICRO}"

perl -pi -e 's/__version__ = \([^"]*\)/__version__ = \($ENV{PYVERSION}\)/' aware/__version__.py
perl -pi -e 's/version = "[^"]*"/version = "$ENV{PYPROJECTTOMLVERSION}"/' pyproject.toml
