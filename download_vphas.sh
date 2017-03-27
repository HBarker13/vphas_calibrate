#! /bin/bash
#Download files from vphas website using the request code.  Rename using vphas pointing
# $1 = reqNo, $2=vphas_num

wget --http-user=vphas --http-passwd=vPhaS4vsT http://casu.ast.cam.ac.uk/~eglez/vstsp/requests/$1/filelist

mv filelist filelist_$2

wget --http-user=vphas --http-passwd=vPhaS4vsT -m -B http://casu.ast.cam.ac.uk/~eglez/vstsp/requests/$1/ -i filelist_$2

echo "Download completed"

mv $PWD/casu.ast.cam.ac.uk/~eglez/vstsp/requests/$1  $PWD/vphas_$2
echo "Downloaded files moved"

rm -r casu.ast.cam.ac.uk
echo "Download directories deleted"

mv filelist_$2 $PWD/vphas_$2/filelist_$2
echo "filelist moved into dir"




