#! /bin/sh
#
# Links all of the unbiased CDM files to the current directory, making
# the level dirs if they don't exist. Then links the CDM velocity to
# the baryon velocity, i.e. same velocity transfer function.
#
# cd to the directory where you want to store the biased IC files
# (e.g. simply "./biased") and set $base to be the path to the
# unbiased directory (from the directory that contains both sets of
# ICs, not the individual IC directory) then run this script to link
# all of the cdm IC files
# 
# So if your directory looked like this:
#
# |-- ics
# |   |-- same_vtf/
# |   |-- unbiased/
#
# you would cd to same_vtf/ and set base to unbiased/

# Make sure you are in the biased directory
base=unbiased
echo "Linking from ../${base}"
read -p "Are you in the directory which will contain level dirs? (y/n) " ans 

fields=("ic_poscx" "ic_poscy" "ic_poscz" "ic_velcx" "ic_velcy" "ic_velcz" "ic_deltab" "ic_refmap" "ic_pvar_00001")
xyz=("x" "y" "z")
levels=$(seq 8 16)

if [[ $ans == "y" ]]; then
    for level in $levels; do
	l=$(( $level + 1000 ))
	src=level_${l:1:3}
	
	# Make the level dirs, if they don't exist
	if [[ ! -d ${src} ]]; then
	    mkdir ${src}
	fi

	cd $src
	
	# Link all the CDM files from $base
	for field in ${fields[@]}; do	    
	    ln -s ../../$base/$src/$field $field
	done

	# Link CDM velocity to baryon velocity
	for i in ${xyz[@]}; do
	    ln -s ic_velc${i} ic_velb${i}
	done

	cd ../
	
    done
fi
