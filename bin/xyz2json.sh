#!/bin/bash

ext=$1
file=$2
directory=$3

## usage: xyz2json.sh <structure_files_extension> sheap_file.txt resfile_directory/ > output_file.txt

dim=$(head -2 $file | tail -1 | awk '{print $2}')

m=2
echo '{'
echo '    "meta": {'
echo '        "name": "SHEAP",'
echo '        "description": "Map generated with SHEAP"'
echo '    },'
echo ''
echo '    "properties": {'
echo '        "sheap1": {'
echo '            "target": "structure",'
echo '            "values": ['$(tail -n+3 $file | awk '{print $'"$m"'}' | sed -e ':a;N;$!ba;s/\n/, /g')']'
echo '        },'
m=$((m+1))
echo '        "sheap2": {'
echo '            "target": "structure",'
echo '            "values": ['$(tail -n+3 $file | awk '{print $'"$m"'}' | sed -e ':a;N;$!ba;s/\n/, /g')']'
echo '        },'
m=$((m+1))
if [ $dim -ge 3 ]
then
echo '        "sheap3": {'
echo '            "target": "structure",'
echo '            "values": ['$(tail -n+3 $file | awk '{print $'"$m"'}' | sed -e ':a;N;$!ba;s/\n/, /g')']'
echo '        },'
m=$((m+1))
fi
m=$((m+4))
echo '        "volume": {'
echo '            "target": "structure",'
echo '            "values": ['$(tail -n+3 $file | awk '{print $'"$m"'}' | sed -e ':a;N;$!ba;s/\n/, /g')']'
echo '        },'
m=$((m+1))
echo '        "energy": {'
echo '            "target": "structure",'
echo '            "values": ['$(tail -n+3 $file | awk '{print $'"$m"'}' | sed -e ':a;N;$!ba;s/\n/, /g')']'
echo '        },'
m=$((m+1))
echo '        "count": {'
echo '            "target": "structure",'
echo '            "values": ['$(tail -n+3 $file | awk '{print $'"$m"'}' | sed -e ':a;N;$!ba;s/\n/, /g')']'
echo '        },'
m=$((m+1))
echo '        "radius": {'
echo '            "target": "structure",'
echo '            "values": ['$(tail -n+3 $file | awk '{print $'"$m"'}' | sed -e ':a;N;$!ba;s/\n/, /g')']'
echo '        }'
echo '    },'
echo ''
echo '    "structures": ['
m=$((m-7))
resfiles=( $(tail -n+3 $file | awk '{print $'"$m"'}' ) )
for f in "${resfiles[@]}"
do
echo '        {'
    resfile=${f%:*:*}
    cabal-sheap $ext json < $directory/$resfile.$ext
    if [ "$f" = "${resfiles[-1]}" ]; then
        echo '        }'
    else
        echo '        },'
    fi
done
echo '    ]'
echo '}'
