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

sed -e "s/__version__ = (\"[0-9]\+\", \"[0-9]\+\", \"[0-9]\+\")/__version__ = \($PYVERSION\)/" -i aware/__version__.py
sed -e "s/version = \"[^\"]*\"/version = \"$PYPROJECTTOMLVERSION\"/" -i pyproject.toml
