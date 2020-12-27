#!/bin/bash                         #-- what is the language of this shell
#                                  #-- Any line that starts with #$ is an instruction to SGE
#$ -S /bin/bash                     #-- the shell for the job
##$ -o job_outputs                        #-- output directory (fill in)
##$ -e job_outputs                        #-- error directory (fill in)
#$ -o /dev/null                        #-- output directory (fill in)
#$ -e /dev/null                        #-- error directory (fill in)
#$ -cwd                            #-- tell the job that it should start in your working directory
#$ -r y                            #-- tell the system that if a job crashes, it should be restarted
#$ -j y                            #-- tell the system that the STDERR and STDOUT should be joined
#$ -l mem_free=5G                  #-- submits on nodes with enough free memory (required)
#$ -l scratch=1G         #-- SGE resources (home and scratch disks)
#$ -l h_rt=24:00:00                #-- runtime limit (see above; this requests 24 hours)
#$ -t 1-25860                        #-- remove first '#' to specify the number of
                                   #-- tasks if desired (see Tips section)

cp -r scaffolds/designable_2lv8_two_LHL_1000_random/ ${TMPDIR}

./scripts/user_scripts_xingjie/match_binding_sites_to_scaffolds.py ${TMPDIR}/designable_2lv8_two_LHL_1000_random/ matching_3res_site_to_2lv8_two_LHL_no_dump/ ${SGE_TASK_ID} -n 3

./scripts/user_scripts_xingjie/match_binding_sites_to_scaffolds.py ${TMPDIR}/designable_2lv8_two_LHL_1000_random/ matching_to_2lv8_two_LHL_no_dump/ ${SGE_TASK_ID} 
