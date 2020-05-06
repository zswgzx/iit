#!/bin/bash

# converts PARREC to nii.
# uses FSL (v 3.3 or 4.0) to do so.
# version 1.8
#  last modified 25 April 2008 to rotate out bvec data so it is patient space.
#  version 1.9.1 April 2009. Modified to include the -Q option
# which includes a rotation matrix in the header, with no rotation.
# so that images are displayed correctly L-R wise in mricron

#version 1.10 May 2009 with modification to read v4.2 output. 
# note this has extra 'number of ASL types' 
# so far have not converted it to read ASL data.

#version 1.12 Oct 2010. Now will (probably) cope when the slices aren't in the
#correct order to start with. Have added something to manually carve up the REC file


##    This program is free software: you can redistribute it and/or modify
##    it under the terms of the GNU General Public License as published by
##    the Free Software Foundation, either version 3 of the License, or
##    (at your option) any later version.

##    This program is distributed in the hope that it will be useful,
##    but WITHOUT ANY WARRANTY; without even the implied warranty of
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##    GNU General Public License for more details.

##    <http://www.gnu.org/licenses/>.


# needed developments
# ?possiblity to have multifile output
# ?use cat instead of avwsplit file may? speed up.
# ?look into using avwcreatehd with -x option to give slope and scale values.
# ? check if no. slices is less than expected from diffusion / dynamics etc
# matrix to cope with different size pixels


Usage() {
    echo ""
    echo "Usage: PARconv input output [-f -q -m]"
    echo "Converts PAR to nifti "
    echo "removes any .PAR .nii .img etc from the filenames you give"
    echo "-f option converts data to 32bit floating point using PAR-REC pixel scaling"
    echo "-q option saves rotation matrix in nifti file (qform matrix)"
    echo "-Q option saves matrix saying file is radiological"
    echo "-m only generates .omat file, no image conversion"
    exit
}
TooMany() {
  echo "Sadly this is a mulitple image file of a type not yet supported"
  echo "bye..."
  export FSLOUTPUTTYPE=$fslot
  exit
}

wrongVersion() {
  echo "This version of PAR REC is not what the tool was designed for."
  echo "Sorry. Bye. Feel free to twiddle at your own risk"
  exit
}

Fine() {
export FSLOUTPUTTYPE=$fslot
exit
}

[ "$1" = "" ] && Usage

# check to see which (if any) version of FSL we have

v3=`fslhd 2>&1 | grep found`
if [ -n "$v3" ]; then
  v2=`avwhd++ 2>&1 | grep found`
  echo $v2
  if [ -n "$v2" ]; then
    echo "can't find a recent version of fsl"
    exit
  fi
  echo v3
  fscreatehd=avwcreatehd++
  fssplit=avwsplit
  fsmerge=avwmerge
  fsmaths=avwmaths++
  fscpgeom=avwcpgeom
  fshd=avwhd++
  fsswapdim=avwswapdim
  fssize=avwsize
else
  echo using fsl v4
  fscreatehd=fslcreatehd
  fssplit=fslsplit
  fsmerge=fslmerge
  fsmaths=fslmaths
  fscpgeom=fslcpgeom
  fshd=fslhd
  fsswapdim=fslswapdim
  fssize=fslsize
fi



tmpfiles=`ls tmp.nii* tmp.img* tmp.hdr* tmp0.nii* tmp0.img* tmp0.hdr*  tmpfile.REC 2>/dev/null`
if [ -n "$tmpfiles"  ]; then  #need " or gets confused with too many arguments
 echo you have files called:
 echo $tmpfiles 
 echo This will confuse things. Please remove
 exit
fi

#set -x

#remove file extensions
#Q1=`${FSLDIR}/bin/remove_ext $1`
#Q1=$1
Q1=`echo $1 | sed s/\.PAR$//`
Q2=`${FSLDIR}/bin/remove_ext $2`


if  [ -a "$Q1".PAR  ]; then 
 echo "$Q1".PAR
else
 echo "$Q1".PAR "is not existing"
 exit
fi
echo $Q2

 imconv=1
 fpt=0
 qf=0
qqf=0
ct=0
for argv in $*; do
  ct=$(($ct+1))
  if [ $ct -gt 2 ]; then
    case $argv in
    "-m") # do matrix only not image conversion
      imconv=0
      ;;
    "-f")
      fpt=1
      ;;
    "-q") 
      qf=1    #output matrix in nifti qform
      ;;
    "-Q")
      qqf=2    #output 'straight' matrix in nifti qform
      ;;
    *)
      echo "warning: unknown option: " $argv 
      ;;
    esac
  fi
done

if [ $qqf -gt 1 ]; then
   if [ $qf -gt 0 ]; then
      echo "-q overrides -Q option - saving angulation matrix"
   else
      qf=2
   fi
fi



if [ $imconv -gt 0 ]; then

if  [ -a "$Q1".REC  ]; then 
 echo "$Q1".REC
else
 echo "$Q1".REC "is not existing"
 exit
fi

tmpfiles=`ls "$Q2".nii* "$Q2".img* "$Q2".hdr 2>/dev/null`
if [ -n "$tmpfiles"  ]; then  #need " or gets confused with too many arguments
 echo "you already have output files:"
 echo $tmpfiles
 echo ". Please remove or rename them"
 exit
fi


#echo  $Q2


  version=`grep "export tool" "$Q1".PAR | awk '{print $8}' | sed s/'\r'//1`
if [ "$version" == "V4.1" ] 
then echo "data version 4.1"
elif [ "$version" == "V4.2" ]
then echo "data version 4.2 : WARNING, STILL EXPERIMENTAL"
else wrongVersion
fi


fslot==$FSLOUTPUTTYPE
export FSLOUTPUTTYPE=NIFTI

 nz=`grep "number of slices" "$Q1".PAR | awk '{print $7}' | sed s/'\r'//1`


stli=`grep -n "sl ec " "$Q1".PAR | sed s/:/\ /g | awk '{print $1+2}' | sed s/'\r'//1`

sinfo=`head -n $stli "$Q1".PAR | tail -1 | sed s/'\r'//1` 



ncp=`grep "number of cardiac" "$Q1".PAR | awk '{print $8}' | sed s/'\r'//1`
nec=`grep "number of echoes" "$Q1".PAR | awk '{print $7}' | sed s/'\r'//1`
ndy=`grep "number of dynamics" "$Q1".PAR | awk '{print $7}' | sed s/'\r'//1`
nmx=`grep "number of mixes" "$Q1".PAR | awk '{print $7}' | sed s/'\r'//1`
nbs=`grep "number of diffusion" "$Q1".PAR | awk '{print $8}' | sed s/'\r'//1`
ngd=`grep "number of gradient" "$Q1".PAR | awk '{print $8}' | sed s/'\r'//1`
tr=`grep "Repetition" "$Q1".PAR | awk '{print $6/1000.0}' | sed s/'\r'//1`
nasl=`grep "Number of label types" "$Q1".PAR | awk '{print $9}' | sed s/'\r'//1`
#echo nasl $nasl
if [ "$version" == "V4.2" ]; then
	[ $nasl -gt 0 ] && TooMany
fi

nim=0
ddim=0
[ $ncp -gt 1 ] && TooMany
#[ $nec -gt 1 ] && TooMany
if [ $ndy -gt 1 ]; then
  [ $(($ncp+$nec+$nmx+$nbs+$ngd)) != 5 ] && TooMany
  nim=$ndy
  ddim=$tr
fi
[ $nmx -gt 1 ] && TooMany
#[ $nbs -gt 1 ] && TooMany
if [ $(($ngd+$nbs)) -gt 2 ]; then
  [ $(($ncp+$nec+$nmx+$ndy)) != 4 ] && TooMany
  nim=$ngd                       #`echo $nbs $ngd | awk '{print ($1-1)*$2+1}'`
  ddim=1
fi
#extract image size orientation from first slice. This is assuming they're all the same...

orient=`echo $sinfo | awk '{print $26}'`
#1 =ax ; 3 = cor ; 2 = sag
bits=`echo $sinfo | awk '{print $8}'`
nx=`echo $sinfo | awk '{print $10}'`
ny=`echo $sinfo | awk '{print $11}'`

ddz=`echo $sinfo | awk '{print $23 + $24}'`
ddx=`echo $sinfo | awk '{print $29}'`
ddy=`echo $sinfo | awk '{print $30}'`

cx=`echo $ddx $nx | awk '{print $1*$2/2.0}'`
cy=`echo $ddy $ny | awk '{print $1*$2/2.0}'`
cz=`echo $ddz $nz | awk '{print $1*$2/2.0}'`


case $bits in 
8)
  bt=2;
  ;;
16)
  bt=4;
  ;;
#32)
#  bt=16;
#  ;;
*) echo $bits "bits not yet supported";
  Fine
  ;;
esac


#sort out multiple image types

#extract slice information from PAR file.
nlines=`wc "$Q1".PAR | awk '{print $1}' | sed s/'\r'//1`
imlines=`echo $nlines $stli | awk '{print $1-$2+1}'`
tail "$Q1".PAR -n $imlines | head -n -2 > tmp.partxt
#check all the slices have same number of bits
  nonmatch1=` awk '{print $8}' tmp.partxt | grep -v -e "$bits" | sed s/'\r'//1`
  nonmatch2=` awk '{print $10}' tmp.partxt | grep -v -e "$nx" | sed s/'\r'//1`
  nonmatch3=` awk '{print $11}' tmp.partxt | grep -v -e "$ny" | sed s/'\r'//1`
  nonmatch=`echo $nonmatch1 $nonmatch2 $nonmatch3`
  if [ -n "$nonmatch" ]; then
       echo "not all the images have same size"
       echo "aborting"
       Fine
  fi
slicesize=`echo $bt $nx $ny | awk '{print $1*$2*$3/2}'`
echo $slicesize

#check for multiple vols
  nvols=`echo $imlines $nz | awk '{print ($1 - 2)/$2}'`
  nint=`echo $imlines $nz | awk '{print ($1 - 2)%$2}'`

if [ $nint -gt 0 ]; then
 echo "total number of lines in PAR file not a multiple of nz"
 echo "exiting"
 Fine
fi



  sorted=`sort -k 43 -k 3 -k 2 -k 1 -n -c tmp.partxt 2>&1 | grep disorder`  # pipe stderr output of sort into grep
#check to see if the file is sorted correctly to start with.

  if [ -n "$sorted" ]; then
    echo "slices and dynamics not sorted. This may be too difficult for me..."
    echo "trying a manual reordering of the REC file..."
#    echo "bye"
#    Fine
   mv tmp.partxt tmp1.partxt
   sort  -k 43 -k 42 -k 3 -k 2 -k 1 -n  tmp1.partxt > tmp2.partxt
#the -k 42 sorts by diffusion number which hopefully puts the extra calculated 
#MD image at the end
## split the file up
   split -a 5 -b $slicesize -d "$Q1".REC  mjfxxx
combine="cat "
## and stick it back together again.
echo sticking together
impos=` awk '{printf " mjfxxx%05d ",$7}' tmp2.partxt`
combine=`echo $combine $impos`

#   for i in `seq 1 $(($imlines-2))`; do
#   iminus=$(($i-1))


  $combine > tmpfile.REC
  rm mjfxxx*
#redo par file with ordering of images in the new REC file.
  awk 'BEGIN {ORS=" ";j=0} { for (k=1;k<=6;k++) print " " $k;print j; for (k=8;k<=NF;k++) print $k; print RS;j++} ' tmp2.partxt >>tmp3.partxt
  tlen=`wc tmp2.partxt | awk '{print $1}'`
  head -n $tlen tmp3.partxt > tmp.partxt
## these seemingly unnecessary lines make the file the same size.
# other wise, tail -1 does not work correctly
  fi  


# if multiple volumes
  if [ "(" $ngd -gt 1 -o $ndy -gt 1 ")" ]; then
    if [ $nvols -gt $nim ]; then
      if [ $ngd -gt 1 ]; then
	if [ $(( $nvols - 1 )) -eq $ngd ]; then 
          # if first 3 gradients are othogonal, last volume is calculated isotropic DWI
          # which we don't really want, so check if this seems to be the case, and ignore
          tmpbvec=`tail -1 tmp.partxt | awk '{printf("%d",100*($48+$46+$47))}' | sed s/'\r'//1`
          tmpbval=`tail -1 tmp.partxt | awk '{printf("%d",$34)}' | sed s/'\r'//1`
          if [ "(" $tmpbvec -eq 0 -a $tmpbval -gt 0 ")" ]; then
            echo "last volume seems to be a calculated isotropic diffusion image. Ignoring this"
#		  if [ -n "$sorted" ]; then
#			echo "DATA WAS ALSO UNSORTED. I HAVE NO IDEA WHAT THIS DOES TO THE ISOTROPIC IMAGE"
#		fi

            # remove details of errant volume
            nvols=$ngd
            head -n -$nz tmp.partxt > tmp0.partxt
            mv tmp0.partxt tmp.partxt
          else
            echo "Strange set of diffusion coeffs."
            Fine
          fi
        else
        echo  "multiple echo/thing and diffusion/dynamics not supported"
        Fine
        fi
      else
      echo
      fi
    else
     nvols=1
    fi
  fi


  prss=`echo $sinfo | awk '{print $14}'`
  prrs=`echo $sinfo | awk '{print $13}'`
  prri=`echo $sinfo | awk '{print $12}'`
  nonmatch1=` awk '{print $14}' tmp.partxt | grep -v -e "$prss" | sed s/'\r'//1`
  nonmatch2=` awk '{print $12}' tmp.partxt | grep -v -e "$prri" | sed s/'\r'//1`
  nonmatch3=` awk '{print $13}' tmp.partxt | grep -v -e "$prrs" | sed s/'\r'//1`
  nomatch=`echo $nonmatch1 $nonmatch2 $nonmatch3`



if [ $nvols -gt 1 ]; then
sort -k 43 -k 2 -k 6 -k 5 -k 3 -k 1 -n  tmp.partxt > tmp0.partxt
nim=$(( $nvols * $nz ))
if [ $nim -gt 9999 ]; then
  echo "you have too much data"
  Fine
fi
ddim=1
$fscreatehd $nx $ny 1 $nim $ddx $ddy $ddz $ddim 0 0 0 $bt tmp0 
$fscreatehd $nx $ny $nz $nvols $ddx $ddy $ddz $ddim 0 0 0 $bt tmpg
#gunzip tmp0.nii.gz


head -c 352 tmp0.nii >  tmp.nii
if [ -n "$sorted" ]; then
  cat tmpfile.REC >> tmp.nii
else	
  cat "$Q1".REC >> tmp.nii
fi
$fssplit tmp.nii
rm tmp.nii*
for i in `seq 1 $(($nvols/$ndy))`;
    do
    h1=$(($i*$nz*$ndy))
    
    ndynz=$(($ndy*$nz))

     head -n $h1 tmp0.partxt | tail -n $ndynz > tmpz.partxt 
     sinfo=`head -n 1 tmpz.partxt | sed s/'\r'//1`
     prss=`echo $sinfo | awk '{print $14}'`
     prrs=`echo $sinfo | awk '{print $13}'`
     prri=`echo $sinfo | awk '{print $12}'`
     nonmatch1=` awk '{print $14}' tmpz.partxt | grep -v -e "$prss" | sed s/'\r'//1`
     nonmatch2=` awk '{print $12}' tmpz.partxt | grep -v -e "$prri" | sed s/'\r'//1`
     nonmatch3=` awk '{print $13}' tmpz.partxt | grep -v -e "$prrs" | sed s/'\r'//1`
     nonmatch=`echo $nonmatch1 $nonmatch2 $nonmatch3`
     if [ -n "$nonmatch" ]; then
       echo "not all the images have same multiplication factor"
       echo "aborting"
       Fine
     fi
     if [ $ngd -gt 1 ]; then
      cat tmp.partxt | awk '{if ($1 == 1) print $48,-1.0*$46,$47}' > tmp_bvecs.txt #"$Q2"_bvecs.txt
      cat tmp.partxt | awk '{if ($1 == 1) print $34}' > "$Q2"_bvals.txt
     fi

     pscale=`echo $sinfo | awk '{print 1.0/$14}'`
     poffs=`echo $sinfo | awk '{print $12/$13/$14}'`

# merge all the slices for this volume into a file
     zfiles=`awk '{printf(" vol%04d" , $7) }' tmpz.partxt `
     $fsmerge -z tmpz"$i" $zfiles

# multiply by scaling if reqd
     if [ $fpt == 1 ]; then
       $fsmaths tmpz"$i" -mul $pscale -add $poffs mtmpz"$i" -odt float
     else
       $fsmaths tmpz"$i" mtmpz"$i" 
     fi
#merge into one multivol file
     if [ $i == 1 ]; then
       $fsmerge -t tmp mtmpz"$i"
     else
       $fsmerge -t tmp tmp mtmpz"$i"
     fi
#     imrm tmpz"$i" mtmpz"$i"
   done

#get geometry right ie nx ny nz nvols
# rather than nx ny 1 nz*nvols
   $fscpgeom tmpg tmp
   
   rm vol*.nii
   rm tmp0.partxt tmpz.partxt
   rm tmpz*.nii* mtmpz*.nii* tmpg.nii*
else
# either single scan or diffusion only or dynamic only

#check that for diffusion and dynamic scans they all have same scaling

  if [ -n "$nomatch" ]; then
    echo "not all the images have same multiplication factor"
    echo "aborting"
    Fine
  fi
  if [ $ngd -gt 1 ]; then
   cat tmp.partxt | awk '{if ($1 == 1) print $48,-1.0*$46,$47}' > tmp_bvecs.txt #"$Q2"_bvecs.txt
   cat tmp.partxt | awk '{if ($1 == 1) print $34}' > "$Q2"_bvals.txt
  fi


# floating point scale and offset. (nb mricro (sept 2007) gives 'displayed value' not 'floating point'
  pscale=`echo $sinfo | awk '{print 1.0/$14}'`
  poffs=`echo $sinfo | awk '{print $12/$13/$14}'`

# generate header to accommodate file
#$fscreatehd $nx $ny $nz $nim $ddx $ddy $ddz $ddim $cx $cy $cz $bt tmp 
  $fscreatehd $nx $ny $nz $nim $ddx $ddy $ddz $ddim 0 0 0 $bt tmp 
#I am setting origin to zero, else you get a rotation matrix with file

#gunzip tmp.nii.gz


  tmphd=`$fshd -x tmp | grep -n "/>"  | sed s/:/\ /g | awk '{print $1-1}'`
  $fshd -x tmp | head -n $tmphd > tmp.hdtxt

# do scaling if required

  if [ $fpt == 1 ]; then
# scale images
    echo scl_slope "$pscale" >> tmp.hdtxt
    echo scl_inter "$poffs"  >> tmp.hdtxt
  fi

  echo "/>" >> tmp.hdtxt

  rm tmp.nii
  $fscreatehd tmp.hdtxt tmp0
#  gunzip tmp0.nii.gz


  head -c 352 tmp0.nii >  tmp.nii
  if [ -n "$sorted" ]; then
    cat tmpfile.REC >> tmp.nii
  else	
    cat "$Q1".REC >> tmp.nii
  fi
  
  rm tmp.hdtxt
fi  #end of multi vs single volume (nvols > 1)

# sort out coronal sagittal so they're in axial format

case $orient in
1) #axial
  $fsswapdim tmp x -y z "$Q2"
  ;;
3) #coronal 
  $fsswapdim tmp x -z -y "$Q2"
  ;;
2) #sagittal
  $fsswapdim tmp -z -x -y "$Q2"
  ;;
*) 
  echo "wierd orientation parameter";
  Fine
  ;;
esac
rm tmp0.nii
rm tmp.nii*

gzip "$Q2".nii

fi # end of image extraction

# make matrix to shift and rotate centre of image to 0,0,0 and no angle.

  AN1ap=`grep "Angulation" "$Q1".PAR | awk '{print $4}' | sed s/'\r'//1`
  AN1fh=`grep "Angulation" "$Q1".PAR | awk '{print $5}' | sed s/'\r'//1`
  AN1rl=`grep "Angulation" "$Q1".PAR | awk '{print $6}' | sed s/'\r'//1`

  TR1ap=`grep "Off Centre" "$Q1".PAR | awk '{print $7}' | sed s/'\r'//1`
  TR1fh=`grep "Off Centre" "$Q1".PAR | awk '{print $8}' | sed s/'\r'//1`
  TR1rl=`grep "Off Centre" "$Q1".PAR | awk '{print $9}' | sed s/'\r'//1`


MD1x=`$fssize "$Q2" | grep dim1 | awk ' BEGIN {sum=1}; {sum *= $2}; END {print sum/2.0}'`
MD1y=`$fssize "$Q2" | grep dim2 | awk ' BEGIN {sum=1}; {sum *= $2}; END {print sum/2.0}'`
MD1z=`$fssize "$Q2" | grep dim3 | awk ' BEGIN {sum=1}; {sum *= $2}; END {print sum/2.0}'`

ddx=`$fssize "$Q2" | grep pixdim1 | awk '{print $2/2}'`
ddy=`$fssize "$Q2" | grep pixdim2 | awk '{print $2/2}'`
ddz=`$fssize "$Q2" | grep pixdim3 | awk '{print $2/2}'`
#echo $MD1x $MD1y $MD1z $ddx $ddy $ddz



#correct for location of pixel middle in translation
#ddx=0
#ddy=0
#ddz=0
MD1x=`echo $MD1x $ddx | awk '{print $1 - $2}'`
MD1y=`echo $MD1y $ddy | awk '{print $1 - $2}'` 
MD1z=`echo $MD1z $ddz | awk '{print $1 - $2}'`

dx=`echo $TR1rl  | awk '{print $1}'`
dy=`echo $TR1ap  | awk '{print -$1}'`
dz=`echo $TR1fh  | awk '{print $1}'`


makerot -t $AN1ap -a 0,-1,0 -c $ddx,$ddy,$ddz > Y1.omat
makerot -t $AN1fh -a 0,0,1 -c  $ddx,$ddy,$ddz > Z1.omat
makerot -t $AN1rl -a 1,0,0 -c  $ddx,$ddy,$ddz > X1.omat
#makerot -t $AN1ap -a 0,-1,0 -c 0,0,0 > Y1.omat
#makerot -t $AN1fh -a 0,0,1 -c  0,0,0 > Z1.omat
#makerot -t $AN1rl -a 1,0,0 -c  0,0,0 > X1.omat

echo 1 0 0 $dx > dxyz.omat
echo 0 1 0 $dy >> dxyz.omat
echo 0 0 1 $dz >> dxyz.omat
echo 0 0 0 1   >> dxyz.omat

echo 1 0 0 -"$MD1x" > dxyz1.omat
echo 0 1 0 -"$MD1y" >> dxyz1.omat
echo 0 0 1 -"$MD1z" >> dxyz1.omat
echo 0 0 0 1   >> dxyz1.omat


convert_xfm -omat tmp0.omat -concat Z1.omat dxyz1.omat
convert_xfm -omat tmp1.omat -concat Y1.omat tmp0.omat
convert_xfm -omat tmp2.omat -concat X1.omat tmp1.omat
convert_xfm -omat "$Q2"_to0.omat -concat dxyz.omat tmp2.omat
convert_xfm -omat  "$Q2"_to0_inv.omat -inverse  "$Q2"_to0.omat


if [ $imconv -gt 0 ]; then
   if [ $ngd -gt 1 ]; then

#if diffusion rotate diffusion bvecs to subject plane.

xra=`head -1  "$Q2"_to0_inv.omat | awk '{print $1}'`
xrb=`head -1  "$Q2"_to0_inv.omat | awk '{print $2}'`
xrc=`head -1  "$Q2"_to0_inv.omat | awk '{print $3}'`
yra=`head -2  "$Q2"_to0_inv.omat  |tail -1| awk '{print $1}'`
yrb=`head -2  "$Q2"_to0_inv.omat  |tail -1| awk '{print $2}'`
yrc=`head -2  "$Q2"_to0_inv.omat  |tail -1| awk '{print $3}'`
zra=`head -3  "$Q2"_to0_inv.omat  |tail -1| awk '{print $1}'`
zrb=`head -3  "$Q2"_to0_inv.omat  |tail -1| awk '{print $2}'`
zrc=`head -3  "$Q2"_to0_inv.omat  |tail -1| awk '{print $3}'`
echo '' >  "$Q2"_bvecs.txt
cat tmp_bvecs.txt | awk '{print $1*'"$xra"'+$2*'"$xrb"'+$3*'"$xrc"', $1*'"$yra"'+$2*'"$yrb"'+$3*'"$yrc"', $1*'"$zra"'+$2*'"$zrb"'+$3*'"$zrc"'}'   >> "$Q2"_bvecs.txt
rm tmp_bvecs.txt
fi
fi

#less -X "$Q2".omat

if [ $qf -gt 0 ]; then
# write rotation matrix into nifti file.

  ddx=`$fssize "$Q2" | grep pixdim1 | awk '{print $2}'`
  ddy=`$fssize "$Q2" | grep pixdim2 | awk '{print $2}'`
  ddz=`$fssize "$Q2" | grep pixdim3 | awk '{print $2}'`
  dimx=`$fssize "$Q2" | grep dim1 | grep -v pixdim | awk '{print $2}'`
  dimy=`$fssize "$Q2" | grep dim2 | grep -v pixdim | awk '{print $2}'`
  dimz=`$fssize "$Q2" | grep dim3 | grep -v pixdim | awk '{print $2}'`
  tmphd=`$fshd -x "$Q2" | grep -n "/>"  | sed s/:/\ /g | awk '{print $1-1}'`
  $fshd -x "$Q2" | head -n $tmphd > tmp.hdtxt
  if [ $qf -eq 1 ]; then
  ## use scanner orient matrix
     qfor=`cat "$Q2"_to0.omat`
  else
  ## generate a matrix 1 0 0 dx/2 ; etc for the sake of having a matrix
     qfor=`echo $dimx $dimy $dimz $ddx $ddy $ddz | awk '{print "1 0 0 ",-1*($1-1)*$4/2.0," 0 1 0 ",-1*($2-1)*$5/2.0," 0 0 1 ",-1*($3-1)*$6/2.0," 0 0 0 1"}'`
  fi
#  echo $qfor
  qforpix=`echo $qfor $ddx $ddy $ddz | awk '{print $1*-$17,$2*-$18,$3*-$19,-1*$4,
                                                   $5*$17,$6*$18,$7*$19,$8,
                                                   $9*$17,$10*$18,$11*$19,$12,
                                                   $13,$14,$15,$16}'`
#  echo $qforpix
  echo   "sform_code = '1'" >> tmp.hdtxt
  echo "sto_xyz_matrix = '" $qforpix "'"  >>tmp.hdtxt
  echo "/>" >> tmp.hdtxt
  $fscreatehd tmp.hdtxt tmp0
  $fscpgeom tmp0 "$Q2"  
  rm tmp0.nii
  rm tmp.hdtxt
fi

rm tmp1.omat tmp2.omat tmp0.omat X1.omat Y1.omat Z1.omat dxyz1.omat dxyz.omat
rm tmp.partxt
  if [ -n "$sorted" ]; then
    rm tmpfile.REC tmp1.partxt tmp2.partxt tmp3.partxt
  fi
if [ $imconv -gt 0 ]; then
  if [ "(" -n "$nomatch" -a $nvols -gt 1 -a $fpt == 0 ")" ]; then
    echo "not all the images in a multi volume have same multiplication factor"
    echo "you might want to use -f option"
  fi
fi
