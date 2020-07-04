#!/bin/bash

usage () {
echo "tableview.sh [-l] [-c <max col width>]
   [ -h ]               help
   [ -l ]               don't pipe through \"less\"
   [ -c ]               truncate each column to this width or less"
}

while getopts "hlc:" OPTION; do
     case $OPTION in
         h) usage; exit;;
         l) NOLESS="true";;
         c) MAXCOLS=$OPTARG;;
     esac
done

shift $((OPTIND-1))

MAXCOLS=${MAXCOLS:-0}

INFILE=${1:-}

BASE1="sed 's/ /√/g' $INFILE | sed 's/\t\t/\tø\t/g' | sed 's/\t\t/\tø\t/g'"
BASE2="| column -t | sed 's/√/ /g' | sed 's/ø/ /g'"

if [ "$NOLESS" != "true" ]; then
    ENDBIT="| less -XS"
else
    ENDBIT=""
fi

if [ "$MAXCOLS" != "0" ]; then
    MIDBIT="| cut -f 1-$MAXCOLS"
fi

CMD="$BASE1 $MIDBIT $BASE2 $ENDBIT"

eval $CMD

