#!/bin/sh
## DESCRIPTION: This "bot" automatically inserts jobs to the TLSMD queue by
##              downloading random PDBs from pdb.org
## LAST UPDATE: 2010-01-05
#for i in `cat list_of_pdbids.txt`;
for i in `mysql -B -N -e 'select id from pdb.remarks where nchains<7 and atoms<8400 and tech like "%X-RAY%" and tlsmd IS NULL order by rand() limit 100;'`;
	do
	./run_tlsmd_bot.py $i >>cron_pdbids_`date +%F`.log 2>>cron_pdbids_`date +%F`.err ;
done
