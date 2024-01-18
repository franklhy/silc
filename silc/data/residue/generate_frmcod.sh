#!/usr/bin/env bash
GAFF2=$AMBERHOME/dat/leap/parm/gaff2.dat
awk 'substr($0,0,12)~/n4/ {if (substr($0,3,1)=="-" && substr($0,6,1)=="-" && substr($0,9,1)=="-") print substr($0,0,55)}' ${GAFF2} > tmp


cat tmp | while IFS= read line
do
	item=`echo ${line} | awk '{print substr($0,0,11)}'`
	newitem=`echo ${item} | sed -e 's/n4/ny/g'`
	if [[ `cat ${GAFF2} | grep "${newitem}" | wc -l` -eq 0 ]];
	then
		echo "${line}" | sed -e 's/n4/ny/g' >> new
	fi
done
echo "Dihedrals for ny, copy from n4 dihedrals in gaff2" > frcmod.silc
echo "DIHEDRAL" >> frcmod.silc
cat new >> frcmod.silc
rm tmp
rm new
