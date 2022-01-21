#! /bin/sh

# cd to the directory where you want to store the biased IC files
# (e.g. simply "./biased") and set $base to be the path to the
# unbiased directory (from the directory that contains both sets of
# ICs, not the individual IC directory) then run this script to link
# all of the cdm IC files

# Make sure you are in the biased directory
base=unbiased
echo "Linking from ${base}"
read -p "Are you in the directory containing level dirs? (y/n) " ans 

fields=("ic_poscx" "ic_poscy" "ic_poscz" "ic_velcx" "ic_velcy" "ic_velcz" "ic_refmap" "ic_pvar_00001")
# fields=("ic_pvar_00001")
# levels=("7")
levels=$(seq 8 15)

# ln -s ../$base/level_008/ level_008

if [[ $ans == "y" ]]; then
    for level in $levels; do
	for field in ${fields[@]}; do
	    l=$(( $level + 1000 ))
	    src=level_${l:1:3}
	    cd $src
	    ln -s ../../$base/$src/$field $field
	    # unlink $field
	    cd ../
	done
    done
fi
