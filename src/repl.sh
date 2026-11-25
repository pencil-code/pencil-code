#!/bin/bash


for file in magnetic/*.f90; do
#echo $file
sed -i -e's/character *( *LEN *= *iomsglen *) *, *intent(out) *:: *iomsg/character(LEN=*), intent(out) :: iomsg/' $file

done

