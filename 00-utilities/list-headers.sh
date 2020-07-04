F=$1
if [ "$F" == "" ]; then
    F="-"
fi

paste <(head -n 1 $F | tr '\t' '\n') <(sed -n '2p' $F | tr '\t' '\n') | cat -n | tableview.sh -l
