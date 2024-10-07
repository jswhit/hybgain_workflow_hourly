for diagfile in  ../../C192_hybcov_hourly_esmda1b/2021082922/diag*ensmean*nc4; do
   echo $diagfile
   python checkobtime.py $diagfile
done
